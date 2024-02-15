process CHECKINPUT {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path samplesheet
    val workflow_name

    output:
    path 'samplesheet.valid.csv', emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 
    """
    case '$workflow_name' in
    'Pre Alignment QC')
        echo $workflow_name detected;
        prealnqc.py \\
            $samplesheet \\
            samplesheet.valid.csv
        ;;
    'DNA Alignment QC')
        dnaalnqc.py \\
            $samplesheet \\
            samplesheet.valid.csv
        ;;
    'DNA Alignment')
        dnaaln.py \\
            $samplesheet \\
            samplesheet.valid.csv
        ;;
    'Germline Variant Call')
        germlinevar.py \\
            $samplesheet \\
            samplesheet.valid.csv
        ;;
    *)
        echo -n "Unknown workflow"
        exit 1
        ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}