process SONG_MANIFEST {
    tag "${analysis_id}"
    label 'process_single'

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]
    container "${ params.song_container ?: 'ghcr.io/overture-stack/song-client' }:${ params.song_container_version ?: '5.0.2' }"
    
    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/song-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/song-client/logs"
    }

    input:
    tuple val(meta), val(analysis_id), path(payload), path(upload)

    output:
    path "out/manifest.txt"       , emit: manifest
    path "versions.yml"           , emit: versions
    tuple val(meta), val(analysis_id), path("out/manifest.txt"), path(upload),    emit: manifest_upload

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def song_url = params.song_url_upload ?: params.song_url
    def accessToken = task.ext.api_upload_token ?: "`cat /tmp/rdpc_secret/secret`"
    def VERSION = params.song_container_version ?: '5.0.2'
    def study_id = "${meta.study_id}"
    """
    export CLIENT_SERVER_URL=${song_url}
    export CLIENT_STUDY_ID=${study_id}
    export CLIENT_ACCESS_TOKEN=${accessToken}
    
    sing manifest -a ${analysis_id} -d . -f ./out/manifest.txt $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: ${VERSION}
    END_VERSIONS
    """
}