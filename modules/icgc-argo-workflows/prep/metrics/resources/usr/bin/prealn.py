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

tool_fieldmap = {
  'fastqc': {},
  'cutadapt': {}
}

fra2pct_fields = ['pct_autosomes_15x', 'pct_autosomes_10x', 'pct_autosomes_30x']
integer_fields = ['r1_with_adapters_total', 'r2_with_adapters_total', 'r_with_adapters_total', 'total_sequences', 'sequences_flagged_as_poor_quality']

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

    # aggregate fastqc and cutadapt metrics into sample level based on multiqc data
    if mqc_stats.get('cutadapt'):
      r_with_adapters_total = 0
      r_processed_total = 0
      r_trimmed_total = 0
      r1_with_adapters_total = 0
      r2_with_adapters_total = 0
      pairs_processed_total = 0
      pairs_trimmed_total = 0
      
      for rg_metrics in mqc_stats['cutadapt']:
        if not rg_metrics.get('pairs_processed'):
          r_with_adapters_total += float(rg_metrics['r_with_adapters'])
          r_processed_total += float(rg_metrics['r_processed'])
          r_trimmed_total += float(rg_metrics['r_processed'])*float(rg_metrics['percent_trimmed'])
        else:
          r1_with_adapters_total += float(rg_metrics['r1_with_adapters'])
          r2_with_adapters_total += float(rg_metrics['r2_with_adapters'])
          pairs_processed_total += float(rg_metrics['pairs_processed'])
          pairs_trimmed_total += float(rg_metrics['pairs_processed'])*float(rg_metrics['percent_trimmed'])

      if r_processed_total > 0:
        mqc_stats['metrics'].update({
          'r_with_adapters_total': r_with_adapters_total,
          'pct_trimmed_total': r_trimmed_total / r_processed_total
        })
      else:
        mqc_stats['metrics'].update({
          'r1_with_adapters_total': r1_with_adapters_total,
          'r2_with_adapters_total': r2_with_adapters_total,
          'pct_trimmed_total': pairs_trimmed_total / pairs_processed_total
        })

    if mqc_stats.get('fastqc'):
      total_sequences = []
      sequences_flagged_as_poor_quality = []
      gc_content = []
      qc_status = {}
      qc_metrics = ['basic_statistics', 'per_base_sequence_quality', 
                    'per_sequence_quality_scores', 'per_base_sequence_content', 'per_sequence_gc_content', 
                    'per_base_n_content', 'sequence_length_distribution', 'sequence_duplication_levels', 
                    'overrepresented_sequences', 'adapter_content']
      for fn in qc_metrics:
        qc_status[fn] = set()
      for rg_metrics in mqc_stats['fastqc']:
        total_sequences.append(float(rg_metrics['Total Sequences']))
        sequences_flagged_as_poor_quality.append(float(rg_metrics['Sequences flagged as poor quality']))
        gc_content.append(float(rg_metrics["%GC"])*float(rg_metrics['Total Sequences']))
        for fn in qc_metrics:
          qc_status[fn].add(rg_metrics[fn])
      
      mqc_stats['metrics'].update(
        {
          'total_sequences': sum(total_sequences),
          'sequences_flagged_as_poor_quality': sum(sequences_flagged_as_poor_quality),
          'pct_gc': sum(gc_content) / sum(total_sequences)
        }
      )
      
      for fn in qc_metrics:
        for status in ['fail', 'warning', 'pass']:          
          if status in qc_status[fn]:
            mqc_stats['metrics'].update(
              {fn: status}
            )
            break

    # convert the fraction to percentage for given fields
    for fn in fra2pct_fields:
      if fn not in mqc_stats['metrics']: continue
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
    """
    Python implementation of tool: prep_metrics.py
    """

    parser = argparse.ArgumentParser(description='Tool: prep_metrics')
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
              'cross_contamination_rate': float(row.get('contamination')),
              'cross_contamination_error': float(row.get('error'))
            })

    mqc_stats_updated = {k: v for k, v in mqc_stats.items() if v}

    with open("%s.prealn.argo_metrics.json" % (args.sampleId), 'w') as f:
      f.write(json.dumps(mqc_stats_updated, indent=2))

if __name__ == "__main__":
    main()

