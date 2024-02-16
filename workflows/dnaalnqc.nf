/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
// WorkflowDnaalnqc.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Check input path parameters to see if they exist
def checkPathParamList = [
  params.input, params.fasta, params.fasta_fai, params.fasta_dict, params.bait_interval, params.target_interval,
  params.germline_resource, params.germline_resource_tbi, params.verifybamid2_ud, params.verifybamid2_mu, params.verifybamid2_bed,
  params.autosome_non_gap
]

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Fails when no target_interval file is provided when target==true
if (params.target && !params.target_interval) {
  error("Please provide target_interval file for target-seq data.")
}

// Initialize all input file channels
fasta       = Channel.fromPath(params.fasta).collect()
fasta_fai   = Channel.fromPath(params.fasta_fai).collect()
fasta_dict  = Channel.fromPath(params.fasta_dict).collect()
target_interval  = params.target_interval   ? Channel.fromPath(params.target_interval).collect() : []
bait_interval    = params.target_interval   ? params.bait_interval  ? Channel.fromPath(params.bait_interval).collect() : Channel.fromPath(params.target_interval).collect() : []
germline_resource      = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect() : Channel.value([])
germline_resource_tbi  = params.germline_resource_tbi  ? Channel.fromPath(params.germline_resource_tbi).collect() : Channel.value([])
verifybamid2_resource  = params.verifybamid2_ud ? Channel.fromPath([params.verifybamid2_ud, params.verifybamid2_mu, params.verifybamid2_bed]).collect() : [[], [], []]
intervals    = params.target ? params.target_interval : params.autosome_non_gap

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK                    } from '../subworkflows/local/input_check'
include { CRAM_QC_MOSDEPTH_SAMTOOLS      } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { CRAM_QC_GATK4_CONTAMINATION as CRAM_QC_CALCONT_PAIR   } from '../subworkflows/local/cram_qc_gatk4_contamination/main'
include { CRAM_QC_GATK4_CONTAMINATION_TUMOUR_ONLY as CRAM_QC_CALCONT_TUMOUR_ONLY   } from '../subworkflows/local/cram_qc_gatk4_contamination_tumour_only/main'
include { CRAM_QC_GATK4_CONTAMINATION_TUMOUR_ONLY as CRAM_QC_CALCONT_NORMAL_ONLY   } from '../subworkflows/local/cram_qc_gatk4_contamination_tumour_only/main'
include { PICARD_COLLECTOXOGMETRICS      } from '../modules/local/picard/collectoxogmetrics/main'
include { BAM_QC_PICARD                  } from '../subworkflows/local/bam_qc_picard/main'
// Build intervals if needed
include { PREPARE_INTERVALS              } from '../subworkflows/local/prepare_intervals/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT ARGO MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_INPUT as STAGE_INPUT_ALN } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { STAGE_INPUT as STAGE_INPUT_QC  } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { PAYLOAD_QCMETRICS              } from '../modules/icgc-argo-workflows/payload/qcmetrics/main'
include { PREP_METRICS                   } from '../modules/icgc-argo-workflows/prep/metrics/main'
include { SONG_SCORE_UPLOAD              } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                        } from '../modules/icgc-argo-workflows/cleanup/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC as MULTIQC_ALL      } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { VERIFYBAMID_VERIFYBAMID2    } from '../modules/nf-core/verifybamid/verifybamid2/main'
include { UNTARFILES                  } from '../modules/nf-core/untarfiles/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DNAALNQC {

    ch_versions = Channel.empty()
    ch_reports = Channel.empty()

    // Stage input files
    STAGE_INPUT_ALN(params.study_id, params.analysis_ids, params.input)
    ch_input_sample = STAGE_INPUT_ALN.out.meta_files
    ch_metadata = STAGE_INPUT_ALN.out.meta_analysis
    ch_versions = ch_versions.mix(STAGE_INPUT_ALN.out.versions)

    // don't enable this feature for now
    // // pull qc metrics from other workflows if qc_analysis_ids are provided
    // if (params.qc_analysis_ids) {
    //   STAGE_INPUT_QC(params.study_id, params.qc_analysis_ids, '')
    //   ch_input_qc_files = STAGE_INPUT_QC.out.meta_files
      
    //   ch_input_qc_files.branch {
    //     duplicate_metrics: it[0].qc_tools.split(',').contains('biobambam2:bammarkduplicates2')
    //   }.set{ch_qc_files}
      
    //   // untar the qc tgz file
    //   UNTARFILES(ch_qc_files.duplicate_metrics)

    //   // Gather QC reports
    //   ch_reports  = ch_reports.mix(UNTARFILES.out.files)
    //   // Gather used softwares versions
    //   ch_versions = ch_versions.mix(UNTARFILES.out.versions)
    // }

  
    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai, intervals, params.no_intervals)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined             = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined

    // MOSDEPTH don't need any intervals for WGS
    intervals_for_mosdepth = 
        intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect()

    intervals            = PREPARE_INTERVALS.out.intervals_bed        // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }

    intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map{ intervals, num_intervals ->
        if ( num_intervals < 1 ) [ [], [], num_intervals ]
        else [ intervals[0], intervals[1], num_intervals ]
    }

    //
    // MODULE: Run LOCAL MODULE PICARD_COLLECTOXOGMETRICS
    //
    PICARD_COLLECTOXOGMETRICS (
      ch_input_sample,
      fasta.map{ it -> [[id:it[0].baseName], it] },               // channel: [ val(meta), fasta ]
      fasta_fai.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_fai ]
      fasta_dict.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_dict ]
      intervals_bed_combined
    )

    // Gather QC reports
    ch_reports  = ch_reports.mix(PICARD_COLLECTOXOGMETRICS.out.metrics)       

    // Gather used softwares versions
    ch_versions = ch_versions.mix(PICARD_COLLECTOXOGMETRICS.out.versions)

    //
    // Module: VERIFYBAMID2
    //
    VERIFYBAMID_VERIFYBAMID2 (
      ch_input_sample,
      verifybamid2_resource,
      [],
      fasta
    )

    // Gather QC reports
    ch_reports  = ch_reports.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm)

    // Gather used softwares versions
    ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)

    //
    // LOCAL SUBWORKFLOW: Run Samtools/stats, Mosdepth
    //
    CRAM_QC_MOSDEPTH_SAMTOOLS (
      ch_input_sample,
      fasta,
      intervals_for_mosdepth
    )

    // Gather QC reports
    ch_reports  = ch_reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)
    
    // Gather used softwares versions
    ch_versions = ch_versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    //
    // LOCAL SUBWORKFLOW: Run GATK4_CALCULATECONTAMINATION, GATK4_GATHERPILEUPSUMMARIES, GATK4_GETPILEUPSUMMARIES
    // Logic to separate normal, tumour samples
    ch_input_sample_branch= ch_input_sample.branch {
      normal: it[0].status == 0
      tumour: it[0].status == 1
    }

    // Normal samples
    cram_normal = ch_input_sample_branch.normal.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ]}

    // Tumour samples
    cram_tumour = ch_input_sample_branch.tumour.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ]}
    
    // Tumour only samples
    // 1. Group together all tumour samples by patient ID [ patient1, [ meta1, meta2 ], [ cram1, cram2], [crai1, crai2 ] ]
    cram_tumour_grouped = cram_tumour.groupTuple()
    
    // 2. Join with normal samples, in each channel there is one key per patient now. 
    // Patients without matched normal end up with: [ patient1, [ meta1, meta2 ], [ cram1, cram2], [crai1, crai2 ], meta3, cram3, crai3 ]
    cram_tumour_joined = cram_tumour_grouped.join(cram_normal, failOnDuplicate: true, remainder: true)
    
    // 3. Filter out entries with last entry null
    // Patients without matched normal end up with: [ patient1, [ meta1, meta2 ], [ cram1, cram2], [crai1, crai2 ], null ]
    cram_tumour_joined_filtered = cram_tumour_joined.filter{ it ->  !(it.last()) }
    
    // 4. Transpose [ patient1, [ meta1, meta2 ], [ cram1, cram2], [crai1, crai2 ], null ] 
    // back to [ patient1, meta1, cram1, crai1, null ] [ patient1, meta2, cram2, crai2, null ]
    // and remove patient ID field & null value for further processing [ meta1, cram1, crai1] [ meta2, cram2, crai2]
    cram_tumor_only = cram_tumour_joined_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

    // Tumour - normal pairs
    // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
    cram_pair = cram_normal.cross(cram_tumour)
    .map { normal, tumour  ->
        def meta = [:]
        meta.id = "${tumour[1].sample}_vs_${normal[1].sample}".toString()
        meta.normal_id  = normal[1].sample
        meta.patient    = normal[0]
        meta.sex        = normal[1].sex
        meta.tumour_id   = tumour[1].sample

        tuple([ meta, [ normal[2], tumour[2] ], [ normal[3], tumour[3] ] ])
    }

    CRAM_QC_CALCONT_PAIR (
      cram_pair,
      fasta.map{ it -> [ [ id:'fasta' ], it ] }, // Remap channel to match module/subworkflow
      fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] }, // Remap channel to match module/subworkflow
      fasta_dict.map{ it -> [ [ id:'fasta_dict' ], it ] },
      germline_resource,
      germline_resource_tbi,
      intervals_and_num_intervals
    )

    CRAM_QC_CALCONT_TUMOUR_ONLY (
      cram_tumor_only,
      fasta.map{ it -> [ [ id:'fasta' ], it ] }, // Remap channel to match module/subworkflow
      fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] }, // Remap channel to match module/subworkflow
      fasta_dict.map{ it -> [ [ id:'fasta_dict' ], it ] },
      germline_resource,
      germline_resource_tbi,
      intervals_and_num_intervals 
    )

    CRAM_QC_CALCONT_NORMAL_ONLY (
      ch_input_sample_branch.normal,
      fasta.map{ it -> [ [ id:'fasta' ], it ] }, // Remap channel to match module/subworkflow
      fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] }, // Remap channel to match module/subworkflow
      fasta_dict.map{ it -> [ [ id:'fasta_dict' ], it ] },
      germline_resource,
      germline_resource_tbi,
      intervals_and_num_intervals 
    )

    // Gather QC reports
    ch_reports  = ch_reports.mix(CRAM_QC_CALCONT_PAIR.out.contamination_table)
    ch_reports  = ch_reports.mix(CRAM_QC_CALCONT_NORMAL_ONLY.out.contamination_table)
    ch_reports  = ch_reports.mix(CRAM_QC_CALCONT_TUMOUR_ONLY.out.contamination_table)
    
    // Gather used softwares versions
    ch_versions = ch_versions.mix(CRAM_QC_CALCONT_PAIR.out.versions)
    ch_versions = ch_versions.mix(CRAM_QC_CALCONT_TUMOUR_ONLY.out.versions)
    ch_versions = ch_versions.mix(CRAM_QC_CALCONT_NORMAL_ONLY.out.versions)

    //
    // LOCAL SUBWORKFLOW: Run BAM_QC_PICARD including PICARD_COLLECTHSMETRICS (targeted, nf-core module), 
    // PICARD_COLLECTWGSMETRICS (WGS, local moduel) and PICARD_COLLECTMULTIPLEMETRICS (nf-core module)
    // When update nf-core modules, ensure the local moduel is not accidentally updated.
    ch_bam_bai_bait_target = params.target ?
        ch_input_sample.combine(bait_interval).combine(target_interval) :
        ch_input_sample.map {meta, cram, crai -> [meta, cram, crai, [], []]}
 
    BAM_QC_PICARD (
      ch_bam_bai_bait_target, // channel: [ val(meta), [bam], [bai], [bait_interval], [target_interval]]
      fasta.map{ it -> [[id:it[0].baseName], it] },               // channel: [ val(meta), fasta ]
      fasta_fai.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_fai ]
      fasta_dict.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_dict ]
      intervals_bed_combined
    )

    // Gather QC reports
    ch_reports  = ch_reports.mix(BAM_QC_PICARD.out.coverage_metrics)
    ch_reports  = ch_reports.mix(BAM_QC_PICARD.out.multiple_metrics)
    
    // Gather used softwares versions
    ch_versions = ch_versions.mix(BAM_QC_PICARD.out.versions)

    //
    // NF-Core MODULE: MultiQC
    //
    ch_multiqc = Channel.empty()
    ch_multiqc = ch_multiqc.mix(ch_reports.collect{meta, report -> report}).ifEmpty([])
  
    MULTIQC_ALL (
        ch_multiqc.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions = ch_versions.mix(MULTIQC_ALL.out.versions)

    // Collect Software Versions
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml'))

    // Group QC files by sampleId
    ch_reports
    .transpose()
    .map { meta, files -> [[id: meta.id, study_id: meta.study_id], files] }
    .groupTuple()
    .set{ ch_meta_reports }

    // Parse the multiqc data & qc files
    PREP_METRICS (ch_meta_reports, MULTIQC_ALL.out.data.collect())

    // Combine channels to determine upload status and payload creation
    // make metadata and files match  
    STAGE_INPUT_ALN.out.meta_analysis.map { meta, metadata -> [[id: meta.sample, study_id: meta.study_id], metadata]}
        .unique().set{ ch_meta_metadata }
  
    ch_meta_metadata.join(ch_meta_reports).join(PREP_METRICS.out.metrics_json)
    .set { ch_metadata_files }

    STAGE_INPUT_ALN.out.upRdpc.combine(ch_metadata_files)
    .map{upRdpc, meta, metadata, files, metrics -> 
    [[id: meta.id, study_id: meta.study_id, upRdpc: upRdpc],
      metadata, files, metrics]}
    .branch{
      upload: it[0].upRdpc
    }.set{ch_metadata_files_status}

    // generate payload
    PAYLOAD_QCMETRICS(
        ch_metadata_files_status, CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect()) 

    SONG_SCORE_UPLOAD(PAYLOAD_QCMETRICS.out.payload_files)
    

    // cleanup the files if specified
    if (params.cleanup) {
      // Gather files to remove   
      ch_files = Channel.empty()
      ch_files = ch_files.mix(STAGE_INPUT_ALN.out.meta_files) 
      ch_files.map{ meta, file1, file2 -> [file1, file2]}
      .set { ch_files_to_remove1 }

      PAYLOAD_QCMETRICS.out.payload_files
      .map {meta, payload, files -> files}
      .unique()
      .set { ch_files_to_remove2 }

      ch_files_to_remove = Channel.empty()
      ch_files_to_remove = ch_files_to_remove.mix(MULTIQC_ALL.out.report)
      ch_files_to_remove = ch_files_to_remove.mix(MULTIQC_ALL.out.data)
      ch_files_to_remove = ch_files_to_remove.mix(ch_files_to_remove1)
      ch_files_to_remove = ch_files_to_remove.mix(ch_files_to_remove2)
      CLEANUP(ch_files_to_remove.unique().collect(), SONG_SCORE_UPLOAD.out.analysis_id)    
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
