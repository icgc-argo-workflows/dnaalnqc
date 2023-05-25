
process PREP_SAMPLE {
  tag "${metadata_json.baseName}"
  label 'process_low'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/icgc-argo/dna-seq-processing-tools:base-docker.0.2.1' }"

  input:  // input, make update as needed
    tuple path(metadata_json), path(input)

  output:  // output, make update as needed
    path "out/*sample_sheet.csv", emit: sample_sheet_csv
    path "versions.yml", emit: versions

  script:
    // add and initialize variables here as needed

    """
    mkdir -p out

    main.py \
      -p ${metadata_json} \
      -s ${input} \
      -n ${task.cpus} \
      -o out
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS

    """
}

