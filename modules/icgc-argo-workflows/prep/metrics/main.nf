process PREP_METRICS {
    tag "$meta.sample"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(qc_files)
    path multiqc // optional

    output:
    tuple val(meta), path('*.argo_metrics.json')   , emit: metrics_json
    tuple val(meta), path('*.metrics.json')   , emit: ga4gh_metrics_json, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def workflow_name = workflow.Manifest.name
    
    
    """
    case '$workflow_name' in
    'Pre Alignment QC')
        echo $workflow_name detected;
        prealn.py \\
            -m $multiqc \\
            -s $meta.sample \\
            -q $qc_files
        ;;
    'DNA Alignment QC')
        dnaaln.py \\
            -m $multiqc \\
            -s $meta.sample \\
            -q $qc_files
        ;;
    'RNA Seq Alignment')
        rnaaln.py \\
            -m $multiqc \\
            -s $meta.sample
        ;;
    'Variant Call QC')
        vcfqc.py \\
            -q $qc_files \\
            -s $meta.sample
        ;;
    *)
        echo -n "Unknown workflow"
        exit 1
        ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
