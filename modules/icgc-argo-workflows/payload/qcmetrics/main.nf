
process PAYLOAD_QCMETRICS {
    tag "$meta.id"
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "bioconda::multiqc=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:  // input, make update as needed
      tuple val(meta), path(metadata_analysis), path(files_to_upload), path(multiqc)
      path pipeline_yml

    output:  // output, make update as needed
      tuple val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
      path "versions.yml", emit: versions

    script:
      // add and initialize variables here as needed
      def arg_pipeline_yml = pipeline_yml ? "-p $pipeline_yml" : ''
      def arg_multiqc = multiqc ? "-m $multiqc" : ''
      """
      main.py \
        -f ${files_to_upload} \
        -a ${metadata_analysis} \
        -w "${workflow.manifest.name}" \
        -s ${workflow.sessionId} \
        -v ${workflow.manifest.version} \
        $arg_pipeline_yml \
        $arg_multiqc

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$(python --version | sed 's/Python //g')
      END_VERSIONS
      """
  }
