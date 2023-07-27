//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION        } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATHERPILEUPSUMMARIES   } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES         } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'


workflow CRAM_QC_GATK4_CONTAMINATION_TUMOUR_ONLY {
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

    // Generate pileup summary tables using getepileupsummaries. 
    GETPILEUPSUMMARIES(input_intervals, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)
    
    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_branch = GETPILEUPSUMMARIES.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    pileup_table_to_merge = pileup_table_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()
    
    // Merge Pileup Summaries
    GATHERPILEUPSUMMARIES(pileup_table_to_merge, dict.map{ meta, dict -> [ dict ] })

    // Mix intervals and no_intervals channels together
    pileup_table = Channel.empty().mix(GATHERPILEUPSUMMARIES.out.table, pileup_table_branch.no_intervals)

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table
    CALCULATECONTAMINATION(pileup_table.map{ meta, table -> [ meta, table, [] ] })
    
    // Gather versions of all tools used
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES.out.versions)

    emit:
    pileup_table // channel: [ meta, table ]
    contamination_table    = CALCULATECONTAMINATION.out.contamination    // channel: [ meta, contamination ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation     // channel: [ meta, segmentation ]
    versions = ch_versions // channel: [ versions.yml ]
}
