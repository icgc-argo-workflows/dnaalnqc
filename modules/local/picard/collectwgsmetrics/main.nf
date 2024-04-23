process PICARD_COLLECTWGSMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    path  intervals

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CollectWgsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    def interval_list = intervals
    def intervallist_cmd = ""
    if (intervals =~ /.(bed|bed.gz)$/){
        interval_list = intervals.toString().replaceAll(/.(bed|bed.gz)$/, ".interval_list")
        intervallist_cmd = "picard -Xmx${avail_mem}M  BedToIntervalList --INPUT ${intervals} --OUTPUT ${interval_list} --SEQUENCE_DICTIONARY ${dict} --TMP_DIR ."
    }

    """
    $intervallist_cmd

    picard \\
        -Xmx${avail_mem}M \\
        CollectWgsMetrics \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.CollectWgsMetrics.coverage_metrics \\
        --REFERENCE_SEQUENCE ${fasta} \\
        --INTERVALS $interval_list


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectWgsMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectWgsMetrics.coverage_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectWgsMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
