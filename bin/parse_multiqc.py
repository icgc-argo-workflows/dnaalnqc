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

file_types = ['fastqc', 'cutadapt', 'samtools_stats', 'picard_wgsmetrics', 'picard_OxoGMetrics', 'picard_QualityYieldMetrics']


def get_mqc_stats(multiqc, sampleId):
    mqc_stats = {}
    for f in sorted(glob(multiqc+'/*.txt')):
      for tool_metrics in file_types:
        if f.endswith(tool_metrics+'.txt'):
          with open(f, 'r') as fn: 
            mqc_stats[tool_metrics] = []
            reader = csv.DictReader(fn, delimiter="\t")
            for row in reader:
              if not sampleId in row.get('Sample'): continue
              mqc_stats[tool_metrics].append(row)
                          
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

