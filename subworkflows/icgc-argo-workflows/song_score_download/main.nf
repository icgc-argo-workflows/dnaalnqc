//
// Run SONG/Score clients to download files
//

include { SONG_GET as songGet } from '../../../modules/icgc-argo-workflows/song/get/main'
include { SCORE_DOWNLOAD as scoreDn } from '../../../modules/icgc-argo-workflows/score/download/main'


// please update workflow code as needed
workflow SONG_SCORE_DOWNLOAD {
  take:  // update as needed
    study_analysis      // channel [study_id, analysis_id]

  main:
    ch_versions = Channel.empty()
   
    songGet(study_analysis)
    ch_versions = ch_versions.mix(songGet.out.versions)

    scoreDn(songGet.out.analysis_json)
    ch_versions = ch_versions.mix(scoreDn.out.versions)

  emit:
    analysis_json = songGet.out.json
    files = scoreDn.out.files
    analysis_files = scoreDn.out.analysis_files  // channel: [analysis, [files]]
    versions = ch_versions                     // channel: [ versions.yml ]
}

