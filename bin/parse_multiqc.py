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

fra2pct_fields = ['pct_autosomes_15x', 'pct_autosomes_10x', 'pct_autosomes_30x']

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

    # convert the fraction to percentage for given fields
    for fn in fra2pct_fields:
      if not mqc_stats['metrics'].get(fn): continue
      new_value = round(mqc_stats['metrics'][fn] * 100, 2)
      mqc_stats['metrics'].update({
        fn: new_value
      })

    # aggregate fastqc and cutadapt metrics into sample level based on multiqc data
    if mqc_stats.get('cutadapt'):
      r1_with_adapters_total = 0
      r2_with_adapters_total = 0
      pairs_processed_total = 0
      pairs_trimmed_total = 0
      for rg_metrics in mqc_stats['cutadapt']:
        r1_with_adapters_total += rg_metrics['r1_with_adapters']
        r2_with_adapters_total += rg_metrics['r2_with_adapters']
        pairs_processed_total += rg_metrics['pairs_processed']
        pairs_trimmed_total += rg_metrics['pairs_processed']*rg_metrics['percent_trimmed']

      mqc_stats['metrics'].update({
        'r1_with_adapters_total': r1_with_adapters_total,
        'r2_with_adapters_total': r2_with_adapters_total,
        'percent_trimmed_total': round(pairs_trimmed_total / pairs_processed_total, 2)
      })

    if mqc_stats.get('fastqc'):
      total_sequences = []
      sequences_flagged_as_poor_quality = []
      percent_gc = []
      qc_status = {}
      qc_metrics = ['basic_statistics', 'per_base_sequence_quality', 'per_tile_sequence_quality', 
                    'per_sequence_quality_scores', 'per_base_sequence_content', 'per_sequence_gc_content', 
                    'per_base_n_content', 'sequence_length_distribution', 'sequence_duplication_levels', 
                    'overrepresented_sequences', 'adapter_content']
      for fn in qc_metrics:
        qc_status[fn] = set()
      for rg_metrics in mqc_stats['fastqc']:
        total_sequences.append(rg_metrics['Total Sequences'])
        sequences_flagged_as_poor_quality.append(rg_metrics['Sequences flagged as poor quality'])
        percent_gc.append(rg_metrics["%GC"])
        for fn in qc_metrics:
          qc_status[fn].add(rg_metrics[fn])
      
      mqc_stats['metrics'].update(
        {
          'total_sequences': sum(total_sequences) / len(total_sequences),
          'sequences_flagged_as_poor_quality': sum(sequences_flagged_as_poor_quality) / len(sequences_flagged_as_poor_quality),
          'percent_gc': sum(percent_gc) / len(percent_gc)
        }
      )
      
      for fn in qc_metrics:
        for status in ['fail', 'warning', 'pass']:          
          if status in qc_status[fn]:
            mqc_stats['metrics'].update(
              {fn: status}
            )
            break

    return mqc_stats


def main():
    """
    Python implementation of tool: payload-gen-qc
    """

    parser = argparse.ArgumentParser(description='Tool: multiqc-parser')
    parser.add_argument("-s", "--sampleId", dest="sampleId", required=True,
                        help="Input sampleId", type=str)
    parser.add_argument("-m", "--multiqc", dest="multiqc", required=True, help="multiqc files folder")
    parser.add_argument("-q", "--qc_files", dest="qc_files", required=True, type=str, nargs="+", help="qc files")

    args = parser.parse_args()
    
    # get tool_specific & aggregated metrics from multiqc
    mqc_stats = {}
    if args.multiqc:
      mqc_stats = get_mqc_stats(args.multiqc, args.sampleId)

    # get tool_specific & aggregated metrics from qc_files when they're not retrieved by multiqc
    for fn in sorted(args.qc_files):
      # GATK4 calculateContamination
      if fn.endswith('contamination.table'):
        file_type = 'gatk_contamination'
        if file_type in mqc_stats: continue
        mqc_stats[file_type] = []
        with open(fn, 'r') as f:     
          reader = csv.DictReader(f, delimiter="\t")
          for row in reader:
            mqc_stats[file_type].append(row)
            mqc_stats['metrics'].update({
              'cross_contamination_rate': row.get('contamination'),
              'cross_contamination_error': row.get('error')
            })

    with open("%s.multiqc_data.json" % (args.sampleId), 'w') as f:
      f.write(json.dumps(mqc_stats, indent=2))



if __name__ == "__main__":
    main()

