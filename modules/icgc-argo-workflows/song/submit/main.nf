
process SONG_SUBMIT {
    tag "$meta.study_id"
    label 'process_single'

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]

    container "${ params.song_container ?: 'ghcr.io/overture-stack/song-client' }:${ params.song_container_version ?: '5.0.2' }"
    
    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/song-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/song-client/logs"
    }

    input:
    tuple val(meta), path(payload), path(files)

    output:
    stdout emit: analysis_id
    tuple val(meta), stdout, path(payload), path(files), emit: analysis_files
    path "versions.yml"           , emit: versions

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

    set -euxo pipefail
    sing submit -f ${payload} $args | jq -er .analysisId | tr -d '\\n'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: ${VERSION}
    END_VERSIONS

    """
}