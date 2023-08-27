process MULTIQC_PARSE {
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
    tuple val(meta), path('*.json')   , emit: multiqc_json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 
    

    """
    parse_multiqc.py \\
        -m $multiqc \\
        -s $meta.id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
