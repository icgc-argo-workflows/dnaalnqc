//
// Check input samplesheet and get read channels
//

//include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    samplesheet
        .splitCsv ( header:true, sep:',' )
        .map { create_input_channel(it) }
        .set { reads }

    index = SAMTOOLS_INDEX (reads).bai
    reads_index = []
    reads_index = reads.join(index)

    emit:
    reads_index                                     // channel: [ val(meta), [ reads, reads_index ] ]
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def create_input_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id    = row.biosample_id
    meta.status = 1
    reads_meta = [meta, file(row.bam_cram)]

    // reads_index = SAMTOOLS_INDEX (reads_meta).bai
    // // add path(s) of the fastq file(s) to the meta map
    // reads_index_meta = []
    // reads_index_meta = reads_meta.join(reads_index)

    return reads_meta
}
