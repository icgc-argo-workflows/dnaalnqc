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
    Edmund
"""

import argparse
import json
from glob import glob
import csv
import sys as sys
import re

ga4gh_metrics = [
    "count_snvs",
    "count_insertions",
    "count_deletions",
    "ratio_insertion_deletion",
    "ratio_heterozygous_homzygous_snv",
    "ratio_heterozygous_homzygous_indel",
    "ratio_transitions_transversions_snv",
]

argo_metric_mapping = {
    "snp":"count_snvs",
    "ins":"count_insertions",
    "del":"count_deletions",
    "het_indel":"count_heterozygous_indels",
    "het_snp":"count_heterozygous_snvs",
    "homo_indel":"count_homozygous_indels",
    "homo_snp":"count_homozygous_snvs",
    "ratio_insertion_deletion":"ratio_insertion_deletion",
    "ratio_heterozygous_homzygous_snv":"ratio_heterozygous_homzygous_snv",
    "ratio_heterozygous_homzygous_indel":"ratio_heterozygous_homzygous_indel",
    "ratio_transitions_transversions_snv":"ratio_transitions_transversions_snv", 
}


def get_qc_stats(qc_file,qc_stats):
    if "metrics" not in qc_stats.keys():
        qc_stats['metrics']={}

    if qc_file.endswith(".vcf"):
        qc_metric=qc_file.split("/")[-1].split(".")[-2]
        with open(qc_file, 'r') as file:
            contents = file.read()
            line_count = contents.count('\n')
        qc_stats['metrics'][argo_metric_mapping[qc_metric]]=line_count
    elif qc_file.endswith("bcftools_stats.txt"):
        with open(qc_file, 'r') as file:
            for line in file:
                ###Return TSTV line
                if bool(re.search("^TSTV", line)):
                    qc_stats['metrics']["ratio_transitions_transversions_snv"]=float(line.split("\t")[4]) # return ts/tv of [id, ts , tv, ts/tv, ts (1st ALT), tv (1st ALT), ts/tv (1st ALT)]
    else:
        print("UNRECOGNIZED FILETYPE %s" % (qc_file))
        sys.exit(1)

    ###Calculate ratios

    if argo_metric_mapping["ins"] in qc_stats['metrics'].keys() and argo_metric_mapping["del"] in qc_stats['metrics'].keys() and "ratio_insertion_deletion" not in qc_stats['metrics'].keys():
        if qc_stats['metrics'][argo_metric_mapping["del"]]==0:
            qc_stats['metrics']["ratio_insertion_deletion"]=None
        else:
            qc_stats['metrics']["ratio_insertion_deletion"]=qc_stats['metrics'][argo_metric_mapping["ins"]]/qc_stats['metrics'][argo_metric_mapping["del"]]

    if argo_metric_mapping["het_snp"] in qc_stats['metrics'].keys() and argo_metric_mapping["homo_snp"] in qc_stats['metrics'].keys() and "ratio_heterozygous_homzygous_snv" not in qc_stats['metrics'].keys():
        if qc_stats['metrics'][argo_metric_mapping["homo_snp"]]==0:
            qc_stats['metrics']["ratio_heterozygous_homzygous_snv"]=None
        else:
            qc_stats['metrics']["ratio_heterozygous_homzygous_snv"]=qc_stats['metrics'][argo_metric_mapping["het_snp"]]/qc_stats['metrics'][argo_metric_mapping["homo_snp"]]

    if argo_metric_mapping["het_indel"] in qc_stats['metrics'].keys() and argo_metric_mapping["homo_indel"] in qc_stats['metrics'].keys() and "ratio_heterozygous_homzygous_indel" not in qc_stats['metrics'].keys():
        if qc_stats['metrics'][argo_metric_mapping["homo_indel"]]==0:
            qc_stats['metrics']["ratio_heterozygous_homzygous_indel"]=None
        else:
            qc_stats['metrics']["ratio_heterozygous_homzygous_indel"]=qc_stats['metrics'][argo_metric_mapping["het_indel"]]/qc_stats['metrics'][argo_metric_mapping["homo_indel"]]

def main():
    parser = argparse.ArgumentParser(description='Tool: prep_metrics')
    parser.add_argument("-s", "--sampleId", dest="sampleId", required=True, help="Input sampleId", type=str)
    parser.add_argument("-q", "--qc_files", dest="qc_files", required=True, type=str, nargs="+", help="qc files")

    args = parser.parse_args()

    # get tool_specific & aggregated metrics from multiqc
    qc_stats = {}
    qc_stats['id']=args.sampleId
    for fn in sorted(args.qc_files):
        print(fn)
        get_qc_stats(fn, qc_stats)
    print(qc_stats)
    
    with open("%s.vcfqc.argo_metrics.json" % (args.sampleId), 'w') as f:
      f.write(json.dumps(qc_stats, indent=2))

    ga4gh_qc_dict = {
      'biosample': {
            'id': args.sampleId
        },
      'vcf_qc_metrics': {}
    }

    for k,v in qc_stats.get('metrics', None).items():
        if not k in ga4gh_metrics: continue
        ga4gh_qc_dict['vcf_qc_metrics'].update({k: v}) 

    if ga4gh_qc_dict['vcf_qc_metrics']:
      with open("%s.vcfqc.metrics.json" % (args.sampleId), 'w') as f:
        f.write(json.dumps(ga4gh_qc_dict, indent=2))

if __name__ == "__main__":
    main()
