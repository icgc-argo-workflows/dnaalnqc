/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                    } from '../subworkflows/local/input_check'
include { STAGE_INPUT as STAGE_INPUT_ALN } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { STAGE_INPUT as STAGE_INPUT_QC  } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { PAYLOAD_QCMETRICS              } from '../modules/icgc-argo-workflows/payload/qcmetrics/main'
include { SONG_SCORE_UPLOAD              } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS      } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { CRAM_QC_GATK4_CONTAMINATION as CRAM_QC_CALCONT_PAIR   } from '../subworkflows/local/cram_qc_gatk4_contamination/main'
include { CRAM_QC_GATK4_CONTAMINATION_TUMOUR_ONLY as CRAM_QC_CALCONT_TUMOUR_ONLY   } from '../subworkflows/local/cram_qc_gatk4_contamination_tumour_only/main'
//include { CRAM_QC_GATK4_CONTAMINATION_TUMOUR_ONLY as CRAM_QC_CALCONT_NORMAL_ONLY   } from '../subworkflows/local/cram_qc_gatk4_contamination_tumour_only/main'
include { PICARD_COLLECTOXOGMETRICS      } from '../modules/local/picard/collectoxogmetrics/main'
include { BAM_QC_PICARD               } from '../subworkflows/local/bam_qc_picard/main'

// Build intervals if needed
include { PREPARE_INTERVALS              } from '../subworkflows/local/prepare_intervals/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC as MULTIQC_ALL      } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_T        } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_N        } from '../modules/nf-core/multiqc/main'
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
    // ch_reports_normal = Channel.empty()
    // ch_reports_tumour = Channel.empty()

    // Read in samplesheet, validate and stage input files
    if ( params.local_mode ) {
      if (params.input) {
        ch_input = Channel.fromPath(params.input)
        ch_input_sample = INPUT_CHECK (ch_input).reads_index
      } 
      else { exit 1, 'Input samplesheet must be specified for local mode!' }
    } else if (params.study_id && params.analysis_id_tumour) {
      ch_input = Channel.of([params.study_id, params.analysis_id_tumour])

      if (params.analysis_id_normal) {
        ch_input = ch_input.concat(Channel.of([params.study_id, params.analysis_id_normal]))
      }
      STAGE_INPUT_ALN(ch_input)
      ch_input_sample = STAGE_INPUT_ALN.out.sample_files
      ch_meta_analysis = STAGE_INPUT_ALN.out.meta_analysis
      ch_versions = ch_versions.mix(STAGE_INPUT_ALN.out.versions)

      // pull qc metrics from other workflows if qc_analysis_ids are provided
      if (params.qc_analysis_ids) {
        // study = Channel.of(params.study_id)
        // params.qc_analysis_ids.split(',').view()
        // qc_analysis = Channel.fromList(params.qc_analysis_ids.split(','))
        // ch_input_qc = study.combine(qc_analysis)
        ch_input_qc = [params.study_id, params.qc_analysis_ids]

        STAGE_INPUT_QC(ch_input_qc)
        ch_input_qc_files = STAGE_INPUT_QC.out.sample_files
        
        ch_input_qc_files.branch {
          duplicate_metrics: it[0].qc_tools.split(',').contains('biobambam2:bammarkduplicates2')
        }.set{ch_qc_files}
        
        // untar the qc tgz file
        UNTARFILES(ch_qc_files.duplicate_metrics)

        // Gather QC reports
        ch_reports  = ch_reports.mix(UNTARFILES.out.files.collect{meta, report -> report})
        // Gather used softwares versions
        ch_versions = ch_versions.mix(UNTARFILES.out.versions)

      }
    } else { exit 1, 'study_id & analysis_id_tumour must be specified for rdpc mode!' }

    // // Initialize all input file channels
    fasta       = Channel.fromPath(params.fasta).collect()
    fasta_fai   = Channel.fromPath(params.fasta_fai).collect()
    fasta_dict  = Channel.fromPath(params.fasta_dict).collect()
    bait_interval    = params.bait_interval   ? Channel.fromPath(params.bait_interval).collect() : []
    target_interval  = params.target_interval ? Channel.fromPath(params.target_interval).collect() : []
    germline_resource      = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect() : Channel.value([])
    germline_resource_tbi  = params.germline_resource_tbi  ? Channel.fromPath(params.germline_resource_tbi).collect() : Channel.value([])
    verifybamid2_resource  = Channel.fromPath([params.verifybamid2_ud, params.verifybamid2_mu, params.verifybamid2_bed]).collect()
    intervals    = params.target ? params.target_interval : params.autosome_non_gap
  
    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai, intervals, params.no_intervals)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined             = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined

    // MOSDEPTH don't need any intervals for WGS
    // intervals_for_mosdepth = params.target ?
    //     intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect() :
    //     Channel.value([ [ id:'null' ], [] ])
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
    // MODULE: Run PICARD_COLLECTOXOGMETRICS
    //
    PICARD_COLLECTOXOGMETRICS (
      ch_input_sample,
      fasta.map{ it -> [[id:it[0].baseName], it] },               // channel: [ val(meta), fasta ]
      fasta_fai.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_fai ]
      fasta_dict.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_dict ]
      intervals_bed_combined
    )

    // Gather QC reports
    PICARD_COLLECTOXOGMETRICS.out.metrics
    .branch {
      normal: it[0].status == 0
      tumour: it[0].status == 1
    }
    .set{oxog_metrics}

    ch_reports  = ch_reports.mix(PICARD_COLLECTOXOGMETRICS.out.metrics).view()
    // ch_reports_normal  = ch_reports_normal.mix(oxog_metrics.normal.collect{meta, report -> report})
    // ch_reports_tumour  = ch_reports_tumour.mix(oxog_metrics.tumour.collect{meta, report -> report})
    
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
    VERIFYBAMID_VERIFYBAMID2.out.self_sm
    .branch {
      normal: it[0].status == 0
      tumour: it[0].status == 1
    }
    .set{verifybamid2_metrics}

    ch_reports  = ch_reports.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm)
    // ch_reports_normal  = ch_reports_normal.mix(verifybamid2_metrics.normal.collect{meta, report -> report})
    // ch_reports_tumour  = ch_reports_tumour.mix(verifybamid2_metrics.tumour.collect{meta, report -> report})


    //
    // LOCAL SUBWORKFLOW: Run Samtools/stats, Mosdepth
    //
    CRAM_QC_MOSDEPTH_SAMTOOLS (
      ch_input_sample,
      fasta,
      intervals_for_mosdepth
    )

    // Gather QC reports
    CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports
    .branch {
      normal: it[0].status == 0
      tumour: it[0].status == 1
    }
    .set {mosdepth_samtools_metrics}

    ch_reports  = ch_reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)
    ch_reports.view()
    // ch_reports_normal  = ch_reports_normal.mix(mosdepth_samtools_metrics.normal.collect{meta, report -> report})
    // ch_reports_tumour  = ch_reports_tumour.mix(mosdepth_samtools_metrics.tumour.collect{meta, report -> report})
    
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
    // 1. Group together all tumour samples by patient ID [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ]
    cram_tumour_grouped = cram_tumour.groupTuple()
    
    // 2. Join with normal samples, in each channel there is one key per patient now. 
    // Patients without matched normal end up with: [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ], null ]
    cram_tumour_joined = cram_tumour_grouped.join(cram_normal, failOnDuplicate: true, remainder: true)
    
    // 3. Filter out entries with last entry null
    cram_tumour_joined_filtered = cram_tumour_joined.filter{ it ->  !(it.last()) }
    
    // 4. Transpose [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ] back to [ patient1, meta1, [ cram1, crai1 ], null ] [ patient1, meta2, [ cram2, crai2 ], null ]
    // and remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ]
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

    // Gather QC reports
    ch_reports  = ch_reports.mix(CRAM_QC_CALCONT_PAIR.out.contamination_table)
    ch_reports  = ch_reports.mix(CRAM_QC_CALCONT_PAIR.out.contamination_table_normal)
    ch_reports  = ch_reports.mix(CRAM_QC_CALCONT_TUMOUR_ONLY.out.contamination_table)
    // ch_reports_tumour  = ch_reports_tumour.mix(CRAM_QC_CALCONT_PAIR.out.contamination_table.collect{meta, report -> report})
    // ch_reports_tumour  = ch_reports_tumour.mix(CRAM_QC_CALCONT_TUMOUR_ONLY.out.contamination_table.collect{meta, report -> report})

    // Gather used softwares versions
    ch_versions = ch_versions.mix(CRAM_QC_CALCONT_PAIR.out.versions)
    ch_versions = ch_versions.mix(CRAM_QC_CALCONT_TUMOUR_ONLY.out.versions)

    //
    // SUBWORKFLOW: Run BAM_QC_PICARD including PICARD_COLLECTHSMETRICS (targeted), PICARD_COLLECTWGSMETRICS (WGS) and PICARD_COLLECTMULTIPLEMETRICS
    // 
    ch_bam_bai_bait_target = params.target ?
        ch_input_sample.combine(bait_interval).combine(target_interval) :
        ch_input_sample.map {meta, cram, crai -> [meta, cram, crai, [], []]}

    // ch_input_sample.combine(bait_interval).combine(target_interval) 
    // .set {ch_bam_bai_bait_target}

 
    BAM_QC_PICARD (
      ch_bam_bai_bait_target, // channel: [ val(meta), [bam], [bai], [bait_interval], [target_interval]]
      fasta.map{ it -> [[id:it[0].baseName], it] },               // channel: [ val(meta), fasta ]
      fasta_fai.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_fai ]
      fasta_dict.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_dict ]
      intervals_bed_combined
    )

    // Gather QC reports
    BAM_QC_PICARD.out.coverage_metrics
    .branch {
      normal: it[0].status == 0
      tumour: it[0].status == 1
    }
    .set {coverage_metrics}

    BAM_QC_PICARD.out.multiple_metrics
    .branch {
      normal: it[0].status == 0
      tumour: it[0].status == 1
    }
    .set {multiple_metrics}

    ch_reports  = ch_reports.mix(BAM_QC_PICARD.out.coverage_metrics)
    ch_reports  = ch_reports.mix(BAM_QC_PICARD.out.multiple_metrics)
    // ch_reports_normal  = ch_reports_normal.mix(coverage_metrics.normal.collect{meta, report -> report})
    // ch_reports_normal  = ch_reports_normal.mix(quality_metrics.normal.collect{meta, report -> report})
    // ch_reports_tumour  = ch_reports_tumour.mix(coverage_metrics.tumour.collect{meta, report -> report})
    // ch_reports_tumour  = ch_reports_tumour.mix(quality_metrics.tumour.collect{meta, report -> report})

    // Gather used softwares versions
    ch_versions = ch_versions.mix(BAM_QC_PICARD.out.versions)

    //
    // MODULE: MultiQC
    //
    ch_multiqc = Channel.empty()
    ch_multiqc = ch_multiqc.mix(ch_reports.collect{meta, report -> report}).ifEmpty([])
    // ch_multiqc_tumour = Channel.empty()
    // ch_multiqc_tumour = ch_multiqc_tumour.mix(ch_reports_tumour.collect().ifEmpty([]))
    // ch_multiqc_normal = Channel.empty()
    // ch_multiqc_normal = ch_multiqc_normal.mix(ch_reports_normal.collect().ifEmpty([]))
    // ch_multiqc_normal.collect()

    MULTIQC_ALL (
        ch_multiqc.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions = ch_versions.mix(MULTIQC_ALL.out.versions)

    // MULTIQC_T (
    //     ch_multiqc_tumour.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // ch_versions = ch_versions.mix(MULTIQC_T.out.versions)

    // MULTIQC_N (
    //     ch_multiqc_normal.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // ch_versions = ch_versions.mix(MULTIQC_N.out.versions)

    // Collect Software Versions
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml'))

    // upload QC files and metadata to song/score
    if (!params.local_mode) {
      //
      // Match the QC files with the metadata info
      //
      // ch_meta_analysis
      // .branch {
      //   normal: it[0].status == 0
      //   tumour: it[0].status == 1
      // }.set {ch_meta_analysis_status}

      // ch_meta_analysis_status.normal
      // .combine(ch_multiqc_normal.collect().concat(MULTIQC_N.out.report, MULTIQC_N.out.data).collect().toList())
      // .set {ch_metadata_upload_normal}

      // ch_meta_analysis_status.tumour
      // .combine(ch_multiqc_tumour.collect().concat(MULTIQC_T.out.report, MULTIQC_T.out.data).collect().toList())
      // .set {ch_metadata_upload_tumour}

      ch_reports
      .map { meta, report -> [ [ id: meta.id ], report ]}
      .groupTuple()
      .set {ch_reports_grouped}

      ch_meta_analysis
      .map { meta, analysis -> [ [ id: meta.id ], analysis ]}

      ch_meta_analysis.join(ch_reports_grouped)
      .set { ch_metadata_upload }

      ch_metadata_upload.view()
      // ch_metadata_upload = ch_metadata_upload.concat(ch_metadata_upload_tumour)

      // // generate payload
      // PAYLOAD_QCMETRICS(
      //   ch_metadata_upload, '', '', CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect(), MULTIQC_T.out.data.collect()
      // ) 

    // //SONG_SCORE_UPLOAD(PAYLOAD_QCMETRICS.out.payload_files)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
