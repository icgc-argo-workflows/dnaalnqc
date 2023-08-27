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
    Linda Xiang
"""

import argparse
import json
from glob import glob
import csv

file_types = {
  'fastqc': {},
  'cutadapt': {},
  'samtools_stats': {
     'pct_reads_mapped': 'reads_mapped_percent',
     'pct_reads_properly_paired': 'reads_properly_paired_percent',
     'mean_insert_size': 'insert_size_average',
     'insert_size_std_deviation': 'insert_size_standard_deviation',
     'total_pf_reads': 'sequences',
     'average_base_quality': 'average_quality',
     'average_read_length': 'average_length',
     'reads_duplicated_percent': 'reads_duplicated_percent',
     'non-primary_alignments': 'non-primary_alignments',
     'pairs_on_different_chromosomes': 'pairs_on_different_chromosomes',
     'mismatch_bases_rate': 'error_rate'
  },
  'picard_wgsmetrics': {
     'mean_autosome_coverage': 'MEAN_COVERAGE',
     'pct_autosomes_15x': 'PCT_15X',
     'mad_autosome_coverage': 'MAD_COVERAGE',
     'median_autosome_coverage': 'MEDIAN_COVERAGE',
     'pct_autosomes_10x': 'PCT_10X',
     'pct_autosomes_30x': 'PCT_30X'
  },
  'picard_OxoGMetrics': {
     'oxidation_q_CCG': 'OXIDATION_Q',
     'oxidation_error_rate_CCG': 'OXIDATION_ERROR_RATE'
  },
  'picard_QualityYieldMetrics': {
     'yield_bp_q30': 'PF_Q30_BASES'
  }
}


def get_mqc_stats(multiqc, sampleId):
    mqc_stats = {
       'metrics': {}
    }
    for f in sorted(glob(multiqc+'/*.txt')):
      for tool_metrics in file_types.keys():
        if f.endswith(tool_metrics+'.txt'):
          with open(f, 'r') as fn: 
            mqc_stats[tool_metrics] = []
            reader = csv.DictReader(fn, delimiter="\t")
            for row in reader:
              if not sampleId in row.get('Sample'): continue
              mqc_stats[tool_metrics].append(row)
              for ftype in file_types.keys():
                if not ftype == tool_metrics: continue
                for f1,f2 in file_types[ftype].items():
                  mqc_stats['metrics'][f1] = row.get(f2)
          
                          
    return mqc_stats


def main():
    """
    Python implementation of tool: payload-gen-qc
    """

    parser = argparse.ArgumentParser(description='Tool: multiqc-parser')
    parser.add_argument("-s", "--sampleId", dest="sampleId", required=True,
                        help="Input sampleId", type=str)
    parser.add_argument("-m", "--multiqc", dest="multiqc", required=True, help="multiqc files folder")

    args = parser.parse_args()
    
    # get tool_specific & aggregated metrics from multiqc
    mqc_stats = {}
    if args.multiqc:
      mqc_stats = get_mqc_stats(args.multiqc, args.sampleId)


    with open("%s.multiqc_data.json" % (args.sampleId), 'w') as f:
      f.write(json.dumps(mqc_stats, indent=2))



if __name__ == "__main__":
    main()

