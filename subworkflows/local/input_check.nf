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
}

def create_input_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id    = row.biosample_id
    meta.status = row.status ? row.status : 0
    meta.patient = row.patient ? row.patient : row.biosample_id
    reads_meta = [meta, file(row.bam_cram)]

    return reads_meta
}
