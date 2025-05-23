#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  Ontario Institute for Cancer Research

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Guanqiao Feng
"""

import argparse
import json
from glob import glob
import csv

tool_fieldmap = { # 'target name' : 'original name'
    'hisat2': {
        'paired_total' : 'paired_total',
        'unpaired_total' : 'unpaired_total'
    },
    'star' : {
        'total_reads' : 'total_reads',
        'avg_input_read_length' : 'avg_input_read_length',
        'pct_uniquely_mapped' : 'uniquely_mapped_percent',
        'avg_mapped_read_length' : 'avg_mapped_read_length',
        'num_splices' : 'num_splices',
        'num_annotated_splices' : 'num_annotated_splices',
        'pct_multimapped' : 'multimapped_percent'
    },
    'picard_RnaSeqMetrics': {
        'correct_strand_reads': 'CORRECT_STRAND_READS',
        'pct_ribosomal_bases': 'PCT_RIBOSOMAL_BASES',
        'pct_coding_bases' : 'PCT_CODING_BASES',
        'pct_utr_bases' : 'PCT_UTR_BASES',
        'pct_intronic_bases' : 'PCT_INTRONIC_BASES',
        'pct_intergenic_bases' : 'PCT_INTERGENIC_BASES',
        'pct_mrna_bases' : 'PCT_MRNA_BASES',
        'pct_usable_bases' : 'PCT_USABLE_BASES',
        'pct_correct_strand_reads' : 'PCT_CORRECT_STRAND_READS',
        'median_cv_coverage' : 'MEDIAN_CV_COVERAGE',
        'median_5prime_to_3prime_bias' : 'MEDIAN_5PRIME_TO_3PRIME_BIAS'
    },
    'samtools_stats': {
        'pct_reads_mapped': 'reads_mapped_percent',
        'pct_reads_properly_paired': 'reads_properly_paired_percent',
        'mean_insert_size': 'insert_size_average',
        'insert_size_std_deviation': 'insert_size_standard_deviation',
        'total_pf_reads': 'sequences',
        'average_base_quality': 'average_quality',
        'average_read_length': 'average_length',
        'pct_reads_duplicated': 'reads_duplicated_percent',
        'non-primary_alignments': 'non-primary_alignments',
        'pairs_on_different_chromosomes': 'pairs_on_different_chromosomes',
        'mismatch_bases_rate': 'error_rate'
    }
}

fra2pct_fields = []
integer_fields = ['paired_total', 'unpaired_total', 'total_reads', 'num_splices', 'num_annotated_splices', 'correct_strand_reads']

def get_mqc_stats(multiqc, sampleId):
    mqc_stats = {
        'sample_id': sampleId,
        'metrics': {}
    }
    for f in sorted(glob(multiqc+'/*.txt')):
        for tool_metrics in tool_fieldmap.keys():
            if f.endswith(tool_metrics+'.txt'):
                with open(f, 'r') as fn: 
                    mqc_stats[tool_metrics] = []
                    reader = csv.DictReader(fn, delimiter="\t")
                    for row in reader:
                        if not sampleId in row.get('Sample'): continue
                        mqc_stats[tool_metrics].append(row)
                        for ftype in tool_fieldmap.keys():
                            if not ftype == tool_metrics: continue
                            for f1,f2 in tool_fieldmap[ftype].items():
                                mqc_stats['metrics'][f1] = float(row.get(f2))

    # convert the fraction to percentage for given fields
    for fn in fra2pct_fields:
        if not mqc_stats['metrics'].get(fn): 
            print(f"Field '{fn}' not found in mqc_stats['metrics'] dictionary")
            continue
        new_value = float(mqc_stats['metrics'][fn]) * 100
        mqc_stats['metrics'].update({
            fn: new_value
        })

    # change type to integer for given fields
    for fn in integer_fields:
        if fn not in mqc_stats['metrics']: continue
        new_value = int(mqc_stats['metrics'][fn])
        mqc_stats['metrics'].update({
            fn: new_value
        })

    return mqc_stats

def main():
    parser = argparse.ArgumentParser(description='Tool: prep_metrics')
    parser.add_argument("-s", "--sampleId", dest="sampleId", required=True, help="Input sampleId", type=str)
    parser.add_argument("-m", "--multiqc", dest="multiqc", required=True, help="multiqc files folder", type=str, nargs="+")

    args = parser.parse_args()

    # get tool_specific & aggregated metrics from multiqc
    mqc_stats = {}
    for fn in sorted(args.multiqc):
        mqc_stats = get_mqc_stats(fn, args.sampleId)

    mqc_stats_updated = {k: v for k, v in mqc_stats.items() if v}

    with open("%s.rnaaln.argo_metrics.json" % (args.sampleId), 'w') as f:
        f.write(json.dumps(mqc_stats_updated, indent=2))

if __name__ == "__main__":
    main()
