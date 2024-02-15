
process SCORE_UPLOAD {
    tag "${analysis_id}"
    label 'process_medium'

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]
    container "${ params.score_container ?: 'ghcr.io/overture-stack/score' }:${ params.score_container_version ?: '5.8.1' }"

    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/score-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/score-client/logs"
    }

    input:
    tuple val(meta), val(analysis_id), path(manifest), path(upload)

    output:
    tuple val(meta), val(analysis_id),    emit: ready_to_publish
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def song_url = params.song_url_upload ?: params.song_url
    def score_url = params.score_url_upload ?: params.score_url
    def transport_parallel = params.transport_parallel ?: task.cpus
    def transport_mem = params.transport_mem ?: "2"
    def accessToken = task.ext.api_upload_token ?: "`cat /tmp/rdpc_secret/secret`"
    def VERSION = params.score_container_version ?: '5.8.1'
    """
    export METADATA_URL=${song_url}
    export STORAGE_URL=${score_url}
    export TRANSPORT_PARALLEL=${transport_parallel}
    export TRANSPORT_MEMORY=${transport_mem}
    export ACCESSTOKEN=${accessToken}
    
    score-client upload --manifest ${manifest} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        score-client: ${VERSION}
    END_VERSIONS
    """
}