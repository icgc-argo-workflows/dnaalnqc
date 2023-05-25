//
// Run SONG/Score clients to upload files
//

include { SONG_SUBMIT as songSub } from '../../../modules/icgc-argo-workflows/song/submit/main' 
include { SONG_MANIFEST as songMan } from '../../../modules/icgc-argo-workflows/song/manifest/main' 
include { SCORE_UPLOAD as scoreUp } from '../../../modules/icgc-argo-workflows/score/upload/main' 
include { SONG_PUBLISH as songPub } from '../../../modules/icgc-argo-workflows/song/publish/main' 


workflow SONG_SCORE_UPLOAD {
    take:
        payload_files  //channel: [meta, payload, files]

    main:
        ch_versions = Channel.empty()
        
        // Create new analysis
        songSub(payload_files)
        ch_versions = ch_versions.mix(songSub.out.versions)

        // Generate file manifest for upload
        songMan(songSub.out.analysis_files)
        ch_versions = ch_versions.mix(songMan.out.versions)

        // Upload to SCORE
        scoreUp(songMan.out.manifest_upload)
        ch_versions = ch_versions.mix(scoreUp.out.versions)

        // Publish the analysis
        songPub(scoreUp.out.ready_to_publish)
        ch_versions = ch_versions.mix(songPub.out.versions)

    emit:
        analysis_id = songPub.out.analysis_id
        versions = ch_versions                     // channel: [ versions.yml ]
}