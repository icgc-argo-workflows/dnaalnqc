process CREATE_INTERVALS_BED {
    tag "$intervals"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    path(intervals)

    output:
    path("*.bed")       , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Consecutive intervals are combined into chunks whose total length stays within
    // params.interval_chunk_size (bp); an interval longer than that gets a chunk of its own.
    // Each chunk is written as a 3-column BED file (0-based start).
    //
    // Supported input formats differ only in field separator and start coordinate base:
    //   .bed           chr<TAB>start<TAB>end         0-based start
    //   .interval_list chr<TAB>start<TAB>end<TAB>..  1-based start, '@' header lines are skipped
    //   anything else  chr:start-end                 1-based start
    def name = intervals.toString().toLowerCase()
    def (fs, start_offset) = name.endsWith("bed")           ? ['\t',    0] :
                             name.endsWith("interval_list") ? ['\t',    1] :
                                                              ['[:-]',  1]
    """
    grep -v '^@' ${intervals} | awk -vFS='${fs}' '{
        start = \$2 - ${start_offset}  # BED start (0-based)
        len = \$3 - start
        if (name == "" || (chunk + len) > ${params.interval_chunk_size}) {
            # start a new chunk
            name = sprintf("%s_%d-%d.bed", \$1, start + 1, \$3)
            chunk = 0
        }
        chunk += len
        printf("%s\\t%d\\t%d\\n", \$1, start, \$3) > name
    }'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
