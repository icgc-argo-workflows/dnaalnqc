
include { SONG_SCORE_DOWNLOAD          } from '../../icgc-argo-workflows/song_score_download/main'
include { PREP_SAMPLE                  } from '../../../modules/icgc-argo-workflows/prep/sample/main'
include { CHECKINPUT                   } from '../../../modules/icgc-argo-workflows/checkinput/main'

workflow STAGE_INPUT {

    take:
    study_id // channel: study_id
    analysis_ids // channel: analysis_ids
    samplesheet  // channel: samplesheet
    
    main:
    ch_versions = Channel.empty()

    //If local_mode is specified do not upload To RDPC
    if (params.local_mode){
      upRdpc_flag=false
    } else {
      //Otherwise only upload to RDPC is API_Token is present
      if (params.api_token || params.api_upload_token){
        upRdpc_flag=true
      } else {
        upRdpc_flag=false
      }
    }

    //Apply appropriate action if API_TOKEN is supplied
    if (params.api_token || params.api_download_token){
      //If IDs are present proceed with download otherwise exit
      if (study_id && analysis_ids){

      Channel.from(analysis_ids.split(","))
      .map{analysis_id -> tuple([study_id,analysis_id])}
      .set{ch_study_analysis}

      SONG_SCORE_DOWNLOAD( ch_study_analysis )
      ch_versions = ch_versions.mix(SONG_SCORE_DOWNLOAD.out.versions)

      PREP_SAMPLE ( SONG_SCORE_DOWNLOAD.out.analysis_files )
      ch_versions = ch_versions.mix(PREP_SAMPLE.out.versions)

      analysis_input = PREP_SAMPLE.out.sample_sheet_csv
      } else {
        exit 1, "Using using API_Token, both a study_id and analysis_ids must be specified."
      }
    } else {
      //If no API_Token, check for local samplesheet
      if (samplesheet){
        CHECKINPUT(file(samplesheet,checkIfExists: true),workflow.Manifest.name)
        ch_versions = ch_versions.mix(CHECKINPUT.out.versions)

        analysis_input = CHECKINPUT.out.csv
      } else {
        exit 1, "When no API_TOKEN is provided, a local samplesheet must be provided."
      }
    }
    //Collect meta,data files and analysis_json
    //Two channels for meta,files and meta,analysis_json will be refined afterwards
    analysis_input
    .collectFile(keepHeader: true, name: 'sample_sheet.csv')
    .splitCsv(header:true)
    .map{ row ->
       if (row.analysis_type == "sequencing_experiment" && row.single_end.toLowerCase() == 'false') {
         tuple([
           analysis_type : row.analysis_type,
           id:"${row.sample}-${row.lane}".toString(), 
           study_id:row.study_id,
           patient:row.patient,
           sex:row.sex,
           status:row.status.toInteger(),
           sample:row.sample, 
           read_group:row.read_group.toString(), 
           data_type:'fastq', 
           numLanes:row.read_group_count,
           experiment:row.experiment,
           single_end : row.single_end.toBoolean()
           ], 
           [file(row.fastq_1), file(row.fastq_2)],
           row.analysis_json
           )
       } else if (row.analysis_type == "sequencing_experiment" && row.single_end.toLowerCase() == 'true') {
         tuple([
           analysis_type : row.analysis_type,
           id:"${row.sample}-${row.lane}".toString(), 
           study_id:row.study_id,
           patient:row.patient,
           sex:row.sex,
           status:row.status.toInteger(),
           sample:row.sample, 
           read_group:row.read_group.toString(), 
           data_type:'fastq', 
           numLanes:row.read_group_count,
           experiment:row.experiment,
           single_end : row.single_end.toBoolean()
           ], 
           [file(row.fastq_1)],
           row.analysis_json
           ) 
      } else if (row.analysis_type == "sequencing_alignment") {
        tuple([
          analysis_type : row.analysis_type,
          id:"${row.sample}".toString(),
          study_id:row.study_id,
          patient:row.patient,
          sample:row.sample,
          sex:row.sex,
          status:row.status.toInteger(),
          genome_build:row.genome_build,
          experiment:row.experiment,
          data_type:'cram'], 
          [file(row.cram), file(row.crai)],
          row.analysis_json
          )
      }
      else if (row.analysis_type == "variant_calling") {
        tuple([
          analysis_type : row.analysis_type,
          id:"${row.sample}".toString(),
          study_id:row.study_id, 
          patient:row.patient,
          sample:row.sample,
          sex:row.sex,
          status:row.status.toInteger(), 
          variantcaller:row.variantcaller, 
          genome_build:row.genome_build,
          experiment:row.experiment,
          data_type:'vcf'],
          [file(row.vcf), file(row.tbi)],
          row.analysis_json
          )
      }
      else if (row.analysis_type == "qc_metrics") {
        tuple([
          analysis_type : row.analysis_type,
          id:"${row.sample}".toString(),
          study_id:row.study_id, 
          patient:row.patient,
          sample:row.sample,
          sex:row.sex,
          status:row.status.toInteger(), 
          qc_tools:row.qc_tools,
          genome_build:row.genome_build,
          experiment:row.experiment,
          data_type:'tgz'],
          [file(row.qc_file)],
          row.analysis_json
          )
      }
    }
    .set { ch_input_sample }

    //We want to still have meta when analysis_json doesn't exist
    ch_input_sample.map{ meta,files,analysis ->
      if (analysis){
        tuple([meta,file(analysis)])
      } else {
        tuple([meta,null])
      }
    }
    .unique{it[1]}
    .set{ ch_meta_analysis }

    //Reorganize files as "sequencing_experiment expected input is tuple while other types are flat"
    ch_input_sample.map{ meta,files,analysis ->
      if (meta.analysis_type == "sequencing_experiment"){
        tuple([meta,files])
      } else if (meta.analysis_type == "sequencing_alignment") {
        tuple([meta,files[0],files[1]])
      } else if (meta.analysis_type == "variant_calling") {
        tuple([meta,files[0],files[1]])
      } else if (meta.analysis_type == "qc_metrics") {
        tuple([meta,files[0]])
      }
    }.set{ch_meta_files}

    emit:
    meta_analysis = ch_meta_analysis // channel: [ val(meta), analysis_json]
    meta_files  = ch_meta_files      // channel: [ val(meta), [ files ] ]
    upRdpc = upRdpc_flag
    
    versions = ch_versions                   // channel: [ versions.yml ]
}