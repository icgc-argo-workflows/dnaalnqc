process PREP_METRICS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(qc_files)
    path multiqc

    output:
    tuple val(meta), path('*.argo_metrics.json')   , emit: metrics_json
    tuple val(meta), path('*.metrics.json')   , emit: ga4gh_metrics_json, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 
    

    """
    main.py \\
        -m $multiqc \\
        -s $meta.id \\
        -q $qc_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
