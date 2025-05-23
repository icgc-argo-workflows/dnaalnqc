
include { SONG_SCORE_DOWNLOAD          } from '../../icgc-argo-workflows/song_score_download/main'
include { PREP_SAMPLE                  } from '../../../modules/icgc-argo-workflows/prep/sample/main'
include { CHECKINPUT                   } from '../../../modules/icgc-argo-workflows/checkinput/main'
include { SAMTOOLS_INDEX as BAM_INDEX  } from '../../../modules/icgc-argo-workflows/samtools/index/main'
include { SAMTOOLS_INDEX as CRAM_INDEX } from '../../../modules/icgc-argo-workflows/samtools/index/main'
include { TABIX_TABIX                  } from '../../../modules/icgc-argo-workflows/tabix/tabix/main'

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
    //Collect meta,data files and analysis_json from new samplesheet.csv and handle approrpiately
    analysis_input
    .collectFile(keepHeader: true, name: 'sample_sheet.csv')
    .splitCsv(header:true)
    .map{ row ->
       if (row.analysis_type == "sequencing_experiment" && row.single_end.toLowerCase() == 'false' && row.experiment == "RNA-Seq") {
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
           single_end : row.single_end.toBoolean(),
           library_strandedness : row.library_strandedness
           ], 
           [file(row.fastq_1,checkIfExists: true), file(row.fastq_2,checkIfExists: true)],
           row.analysis_json
           )
       } else if (row.analysis_type == "sequencing_experiment" && row.single_end.toLowerCase() == 'true' && row.experiment == "RNA-Seq") {
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
           single_end : row.single_end.toBoolean(),
           library_strandedness : row.library_strandedness
           ], 
           [file(row.fastq_1,checkIfExists: true)],
           row.analysis_json
           )
       } else if (row.analysis_type == "sequencing_experiment" && row.single_end.toLowerCase() == 'false') {
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
           [file(row.fastq_1,checkIfExists: true), file(row.fastq_2,checkIfExists: true)],
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
           [file(row.fastq_1,checkIfExists: true)],
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
          data_type: "${row.bam_cram}".replaceAll(/^.*\./,"").toLowerCase()], 
          row.bai_crai ? [file(row.bam_cram,checkIfExists: true),file(row.bai_crai,checkIfExists: true)] : [file(row.bam_cram,checkIfExists: true)],
          row.analysis_json
          )
      }
      else if (row.analysis_type == "variant_calling" ) {
        tuple([
          analysis_type : row.analysis_type,
          id:"${row.sample}".toString(),
          study_id:row.study_id, 
          patient:row.patient,
          sample:row.sample,
          sex:row.sex,
          variantcaller:row.variantcaller, 
          genome_build:row.genome_build,
          experiment:row.experiment,
          data_type:'vcf'],
          row.vcf_index ? [file(row.vcf,checkIfExists: true),file(row.vcf_index,checkIfExists: true)] : [file(row.vcf,checkIfExists: true)],
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
          [file(row.qc_file,checkIfExists: true)],
          row.analysis_json
          )
      }
    }
    .set {ch_input_sample}
    
    //Reorganize files as flat tuple except "sequencing_experiment
    ch_input_sample.map{ meta,files,analysis ->
      if (meta.analysis_type == "sequencing_experiment"){
        tuple([meta,files]) //tuple([meta,[read1,read2]])
      } else if (meta.analysis_type == "sequencing_alignment") {
        tuple([meta,files[0],files[1]])
      } else if (meta.analysis_type == "variant_calling") {
        tuple([meta,files])
      } else if (meta.analysis_type == "qc_metrics") {
        tuple([meta,files[0]])
      }
    }.branch{ //identify files that require indexing
      bam_to_index : it[0].analysis_type=='sequencing_alignment' && it[1].size()!=2 && it[0].data_type=='bam'
        return tuple([it[0],it[1]])
      cram_to_index : it[0].analysis_type=='sequencing_alignment' && it[1].size()!=2 && it[0].data_type=='cram'
        return tuple([it[0],it[1]])
      vcf_to_index : it[0].analysis_type=='variant_calling' && it[1].size()!=2
        return tuple([it[0],it[1]])
      indexed : (it[0].analysis_type=='sequencing_alignment' && it[1].size()==2) | (it[0].analysis_type=='variant_calling' && it[1].size()==2)
        return tuple([it[0],it[1][0],it[1][1]])     
      others: (it[0].analysis_type=='sequencing_experiment') | (it[0].analysis_type=='qc_metrics')
        return tuple([it[0],it[1]])
    }.set{ch_index_split}

    //Perform indexiing
    BAM_INDEX(ch_index_split.bam_to_index)
    CRAM_INDEX(ch_index_split.cram_to_index)
    TABIX_TABIX(ch_index_split.vcf_to_index)


    //Combine BAM and BAI into single channel
    ch_index_split.bam_to_index.join(BAM_INDEX.out.bai) //[meta,bam,bai]
    .set{indexed_bam}

    //Combine CRAM and CRAI into single channel
    ch_index_split.cram_to_index.join(CRAM_INDEX.out.crai) //[meta,cram,crai]
    .set{indexed_cram}

    //Combine VCF and TBI into single channel
    ch_index_split.vcf_to_index.join(TABIX_TABIX.out.tbi) //[meta,vcf,tbi]
    .set{indexed_vcf}

    //Combine newly indexed files, previously indexed and others into single channel
    Channel.empty()
    .mix(indexed_bam)
    .mix(indexed_cram)
    .mix(indexed_vcf)
    .mix(ch_index_split.indexed)
    .mix(ch_index_split.others)
    .set{ch_meta_files}


    //We want to still have meta when analysis_json doesn't exist
    ch_input_sample.map{ meta,files,analysis ->
      if (analysis){
        tuple([meta,file(analysis,checkIfExists: true)])
      } else {
        tuple([meta,null])
      }
    }
    .unique{it[1]}
    .set{ ch_meta_analysis }

    emit:
    meta_analysis = ch_meta_analysis // channel: [ val(meta), analysis_json]
    meta_files  = ch_meta_files      // channel: [ val(meta), [ files ] ]
    upRdpc = upRdpc_flag // [boolean]
    
    versions = ch_versions                   // channel: [ versions.yml ]
}