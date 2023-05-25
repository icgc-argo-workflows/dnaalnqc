
include { SONG_SCORE_DOWNLOAD          } from '../../icgc-argo-workflows/song_score_download/main'
include { PREP_SAMPLE                  } from '../../../modules/icgc-argo-workflows/prep/sample/main.nf'

workflow STAGE_INPUT {

    take:
    study_analysis  // channel: study_id, analysis_id

    main:
    ch_versions = Channel.empty()

    SONG_SCORE_DOWNLOAD( study_analysis )
    ch_versions = ch_versions.mix(SONG_SCORE_DOWNLOAD.out.versions)

    PREP_SAMPLE ( SONG_SCORE_DOWNLOAD.out.analysis_files )
    ch_versions = ch_versions.mix(PREP_SAMPLE.out.versions)

    PREP_SAMPLE.out.sample_sheet_csv
    .collectFile(keepHeader: true, name: 'sample_sheet.csv')
    .splitCsv(header:true)
    .map{ row ->
      if (row.analysis_type == "sequencing_experiment") {
        tuple([
          id:"${row.sample}-${row.lane}".toString(), 
          study_id:row.study_id,
          patient:row.patient,
          sex:row.sex,
          status:row.status,
          sample:row.sample, 
          read_group:row.read_group.toString(), 
          data_type:'fastq', 
          size:1, 
          numLanes:row.read_group_count], 
          [file(row.fastq_1), file(row.fastq_2)]) 
      }
      else if (row.analysis_type == "sequencing_alignment") {
        tuple([
          id:"${row.sample}".toString(),
          study_id:row.study_id,
          patient:row.patient,
          sample:row.sample,
          sex:row.sex,
          status:row.status, 
          data_type:'cram'], 
          file(row.cram), file(row.crai))
      }
      else if (row.analysis_type == "variant_calling") {
        tuple([
          id:"${row.sample}".toString(),
          study_id:row.study_id, 
          patient:row.patient,
          sample:row.sample, 
          variantcaller:row.variantcaller, 
          data_type:'vcf'], file(row.vcf), file(row.tbi))
      }
      else if (row.analysis_type == "qc_metrics") {
        tuple([
          id:"${row.sample}".toString(),
          study_id:row.study_id, 
          patient:row.patient,
          sample:row.sample, 
          qc_tools:row.qc_tools, 
          data_type:'tgz'], file(row.qc_file))
      }
    }
    .set { ch_input_sample }

    // ch_input_sample.combine(SONG_SCORE_DOWNLOAD.out.analysis_json)
    // .map { meta, files, analysis_json -> 
    // [meta, analysis_json]
    // }
    // .set { ch_metadata }

    emit:
    analysis_json = SONG_SCORE_DOWNLOAD.out.analysis_json  // channel: [ metadata ] 
    sample_files  = ch_input_sample          // channel: [ val(meta), [ files ] ]
    input_files = SONG_SCORE_DOWNLOAD.out.files // channel: [files]
    versions = ch_versions                   // channel: [ versions.yml ]
}