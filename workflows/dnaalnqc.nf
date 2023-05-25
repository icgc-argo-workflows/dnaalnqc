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
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { STAGE_INPUT as STAGE_INPUT_ALN } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { STAGE_INPUT as STAGE_INPUT_QC } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { PAYLOAD_QCMETRICS } from '../modules/icgc-argo-workflows/payload/qcmetrics/main'
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS   } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { PICARD_COLLECTOXOGMETRICS     } from '../modules/local/picard/collectoxogmetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { BAM_QC_PICARD               } from '../subworkflows/nf-core/bam_qc_picard/main'
include { VERIFYBAMID_VERIFYBAMID2 } from '../modules/nf-core/verifybamid/verifybamid2/main'
include { UNTARFILES } from '../modules/nf-core/untarfiles/main'

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

    // Read in samplesheet, validate and stage input files
    if ( params.local_mode ) {
      if (params.input) {
        ch_input = file(params.input, checkIfExists: true)
        ch_input_sample = INPUT_CHECK (ch_input).reads
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
      } 
      else { exit 1, 'Input samplesheet must be specified for local mode!' }
    } else if (params.study_id && params.analysis_id) {
      ch_input = [params.study_id, params.analysis_id]
      STAGE_INPUT_ALN(ch_input)
      ch_input_sample = STAGE_INPUT_ALN.out.sample_files
      analysis_json = STAGE_INPUT_ALN.out.analysis_json
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
    } else { exit 1, 'study_id & analysis_id must be specified for rdpc mode!' }




    // Initialize all input file channels
    fasta       = Channel.fromPath(params.fasta).collect()
    fasta_fai   = Channel.fromPath(params.fasta_fai).collect()
    fasta_dict  = Channel.fromPath(params.fasta_dict).collect()
    bait_interval    = params.bait_interval   ? Channel.fromPath(params.bait_interval).collect() : []
    target_interval  = params.target_interval ? Channel.fromPath(params.target_interval).collect() : []
    intervals_for_processing = [] 

    //
    // MODULE: Run PICARD_COLLECTOXOGMETRICS
    //
    PICARD_COLLECTOXOGMETRICS (
      ch_input_sample,
      fasta.map{ it -> [[id:it[0].baseName], it] },               // channel: [ val(meta), fasta ]
      fasta_fai.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_fai ]
      intervals_for_processing
    )
    // Gather QC reports
    ch_reports  = ch_reports.mix(PICARD_COLLECTOXOGMETRICS.out.metrics.collect{meta, report -> report})
    // Gather used softwares versions
    ch_versions = ch_versions.mix(PICARD_COLLECTOXOGMETRICS.out.versions)

    //
    // SUBWORKFLOW: Run Samtools/stats, Mosdepth
    //
    CRAM_QC_MOSDEPTH_SAMTOOLS (
      ch_input_sample,
      fasta,
      fasta_fai,
      intervals_for_processing
    )
    // Gather QC reports
    ch_reports  = ch_reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.qc.collect{meta, report -> report})
    // Gather used softwares versions
    ch_versions = ch_versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)


    //
    // SUBWORKFLOW: Run BAM_QC_PICARD including PICARD_COLLECTHSMETRICS (targeted), PICARD_COLLECTWGSMETRICS (WGS)
    // 
    // ch_input_sample.map {meta, cram, crai -> meta}.set {ch_meta}
    // ch_meta.combine(bait_interval).set {ch_bait_interval}
    ch_input_sample.map {meta, cram, crai -> 
      [meta, cram, crai, [], []]
    }
    .set {ch_bam_bai_bait_target}
 
    BAM_QC_PICARD (
      ch_bam_bai_bait_target, // channel: [ val(meta), [bam], [bai], [bait_interval], [target_interval]]
      fasta.map{ it -> [[id:it[0].baseName], it] },               // channel: [ val(meta), fasta ]
      fasta_fai.map{ it -> [[id:it[0].baseName], it] },           // channel: [ val(meta), fasta_fai ]
      fasta_dict.map{ it -> [[id:it[0].baseName], it] }           // channel: [ val(meta), fasta_dict ]
    )
    // Gather QC reports
    ch_reports  = ch_reports.mix(BAM_QC_PICARD.out.coverage_metrics.collect{meta, report -> report})
    ch_reports  = ch_reports.mix(BAM_QC_PICARD.out.multiple_metrics.collect{meta, report -> report})
    // Gather used softwares versions
    ch_versions = ch_versions.mix(BAM_QC_PICARD.out.versions)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    // Match the QC files with the metadata info
    ch_input_sample
    .map { meta, cram, crai -> meta }
    .set{ ch_meta }

    ch_meta.combine(ch_multiqc_files.collect().concat(MULTIQC.out.report, MULTIQC.out.data).collect().toList())
    .combine(analysis_json)
    .set {ch_metadata_upload}


    // Collect Software Versions
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml'))    

    // upload QC files and metadata to song/score
    if (!params.local_mode) {
      // generate payload
      PAYLOAD_QCMETRICS(
        ch_metadata_upload, '', '', CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect()) 

      //SONG_SCORE_UPLOAD(PAYLOAD_QCMETRICS.out.payload_files)
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
