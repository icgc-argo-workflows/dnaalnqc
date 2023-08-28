//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION        } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
// include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION_NORMAL } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATHERPILEUPSUMMARIES_NORMAL   } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATHERPILEUPSUMMARIES_TUMOUR   } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL         } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOUR         } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'


workflow CRAM_QC_GATK4_CONTAMINATION {
    take:
    input                     // channel: [ meta, [input] , [input_index ]]
    fasta                     // channel: [ meta, /path/to/reference/fasta]
    fai                       // channel: [ meta, /path/to/reference/fasta/index]
    dict                      // channel: [ meta, /path/to/reference/fasta/dictionary]
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    intervals                 // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    ch_versions = Channel.empty()

    germline_resource_pileup     = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
    // Move num_intervals to meta map and reorganize channel for MUTECT2_PAIRED module
    .map{ meta, input_list, input_index_list, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input_list, input_index_list, intervals ] }


    pileup = input_intervals.multiMap{  meta, input_list, input_index_list, intervals ->
        tumour: [ meta, input_list[1], input_index_list[1], intervals ]
        normal: [ meta, input_list[0], input_index_list[0], intervals ]
    }

    pileup_normal = pileup.normal.map{ meta, cram, crai, intervals -> [ meta + [ id:meta.normal_id ], cram, crai, intervals ] }
    pileup_tumour = pileup.tumour.map{ meta, cram, crai, intervals -> [ meta + [ id:meta.tumour_id ], cram, crai, intervals ] }

    // Generate pileup summary tables using getepileupsummaries. tumour sample should always be passed in as the first input and input list entries of vcf_to_filter,
    GETPILEUPSUMMARIES_NORMAL(pileup_normal, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)
    GETPILEUPSUMMARIES_TUMOUR(pileup_tumour, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_normal_branch = GETPILEUPSUMMARIES_NORMAL.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_tumour_branch = GETPILEUPSUMMARIES_TUMOUR.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    pileup_table_normal_to_merge = pileup_table_normal_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()
    pileup_table_tumour_to_merge = pileup_table_tumour_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()

    // Merge Pileup Summaries
    GATHERPILEUPSUMMARIES_NORMAL(pileup_table_normal_to_merge, dict.map{ meta, dict ->  dict  })
    GATHERPILEUPSUMMARIES_TUMOUR(pileup_table_tumour_to_merge, dict.map{ meta, dict ->  dict  })

    // remove no longer necessary field: normal_id, tumour_id, num_intervals
    pileup_table_normal = Channel.empty().mix(GATHERPILEUPSUMMARIES_NORMAL.out.table, pileup_table_normal_branch.no_intervals)
        .map{ meta, table -> [ meta - meta.subMap('num_intervals') + [ id:meta.tumour_id + "_vs_" + meta.normal_id ], table ] }

    // remove no longer necessary field: normal_id, tumour_id, num_intervals
    pileup_table_tumour = Channel.empty().mix(GATHERPILEUPSUMMARIES_TUMOUR.out.table, pileup_table_tumour_branch.no_intervals)
        .map{ meta, table -> [ meta - meta.subMap('num_intervals') + [ id:meta.tumour_id + "_vs_" + meta.normal_id ], table ] }

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table
    CALCULATECONTAMINATION(pileup_table_tumour.join(pileup_table_normal, failOnDuplicate: true, failOnMismatch: true))

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table
    // CALCULATECONTAMINATION_NORMAL(pileup_table_normal.map{ meta, table -> [ meta + [id:meta.normal_id], table, [] ] })

    // Gather versions of all tools used

    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_TUMOUR.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_NORMAL.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_TUMOUR.out.versions)

    emit:
    pileup_table_normal // channel: [ meta, table_normal ]
    pileup_table_tumour  // channel: [ meta, table_tumour ]
    contamination_table    = CALCULATECONTAMINATION.out.contamination.map {meta, contamination -> [[id:meta.tumour_id], contamination]}    // channel: [ meta, contamination ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation.map {meta, segmentation -> [[id:meta.tumour_id], segmentation]}     // channel: [ meta, segmentation ]
    // contamination_table_normal    = CALCULATECONTAMINATION_NORMAL.out.contamination.map {meta, contamination -> [[id:meta.normal_id], contamination]}    // channel: [ meta, contamination ]
    // segmentation_table_normal     = CALCULATECONTAMINATION_NORMAL.out.segmentation.map {meta, segmentation -> [[id:meta.normal_id], segmentation]}     // channel: [ meta, segmentation ]
    versions = ch_versions // channel: [ versions.yml ]
}
