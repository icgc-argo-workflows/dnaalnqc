#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import json
import re
import hashlib
from argparse import ArgumentParser
from multiprocessing import cpu_count
import glob
import shutil
import csv

def run_cmd(cmd):
    proc = subprocess.Popen(
                cmd,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
    stdout, stderr = proc.communicate()

    return (
        stdout.decode("utf-8").strip(),
        stderr.decode("utf-8").strip(),
        proc.returncode
    )

def group_readgroup_by_filepair(song_analysis):
    filepair_map_to_readgroup = {}
    # we assume information in song_analysis has gone through
    # all validation checks in sequencing experiment submission
    # note: since advanced SONG validation is not ready, here we still validate uniqueness of
    #       read_group_id_in_bam and submitter_read_group_id
    read_group_id_in_bam_set = set()
    submitter_read_group_id_set = set()
    for rg in song_analysis.get('read_groups'):
        rg['experiment'] = song_analysis['experiment']  # let read group carry experiment
        rg['submitter_sample_id'] = song_analysis['samples'][0]['submitterSampleId']  # let read group carry submitter_sample_id

        file_r1_r2 = (rg.get('file_r1'), rg.get('file_r2'))  # tuple

        if file_r1_r2 not in filepair_map_to_readgroup:
            filepair_map_to_readgroup[file_r1_r2] = {
                'format': 'BAM' if file_r1_r2[0].endswith('.bam') else 'FASTQ',
                'read_groups': [rg]
            }

        else:
            filepair_map_to_readgroup[file_r1_r2]['read_groups'].append(rg)

        # make sure no duplicate of read_group_id_in_bam (when populated) within the same bam
        if rg.get('read_group_id_in_bam'):
            bam_and_rg_id = '%s %s' % (rg.get('file_r1'), rg.get('read_group_id_in_bam'))
            if bam_and_rg_id in read_group_id_in_bam_set:
                sys.exit('Error found: read_group_id_in_bam duplicated in the same BAM: %s' % bam_and_rg_id)
            else:
                read_group_id_in_bam_set.add(bam_and_rg_id)

        # make sure no duplicate of submitter_read_group_id
        if rg['submitter_read_group_id'] in submitter_read_group_id_set:
            sys.exit('Error found: submitter_read_group_id duplicated: %s' % rg['submitter_read_group_id'])
        else:
            submitter_read_group_id_set.add(rg['submitter_read_group_id'])

    return filepair_map_to_readgroup

def filename_to_file(filenames: tuple, files: list) -> tuple:
    name_to_file_map = {}
    for f in files:
      name_to_file_map[os.path.basename(f)] = f

    return (name_to_file_map[filenames[0]], name_to_file_map[filenames[1]] if filenames[1] else None)

def readgroup_id_to_fname(rg_id, input_bam_name='', study_id=None, donor_id=None, sample_id=None):
    friendly_rgid = "".join([ c if re.match(r"[a-zA-Z0-9\.\-_]", c) else "_" for c in rg_id ])
    # use original bam file name and rg_id to calculate the md5sum to avoid new lane bam file name collision
    # use white space (' ') to separate bam name and rg_id
    md5sum = hashlib.md5(("%s %s" % (input_bam_name, rg_id)).encode('utf-8')).hexdigest()

    if not sample_id or not donor_id or not study_id:
        sys.exit('Error: missing study/donor/sample ID in the provided metadata')

    return ".".join([study_id, donor_id, sample_id, friendly_rgid, md5sum])


def generate_fastqs_from_bam(bam, readgroups, cpu=None, sample_sheet=dict(), study_id=None, donor_id=None, sample_id=None, out_dir=None, specimen_id=None, specimen_type=None, tumour_normal_designation=None):
    # get input bam basename, remove extention to use as output subfolder name
    bam_base = os.path.splitext(os.path.basename(bam))[0]
    out_format = bam_base+'/%!.bam'
    if os.path.exists(bam_base) and os.path.isdir(bam_base):
      shutil.rmtree(bam_base)
    os.mkdir(bam_base)
    cmd = ['samtools', 'split', '-f', '%s' % out_format, '-@ %s' % str(cpu), bam]

    stdout, stderr, returncode = run_cmd(" ".join(cmd))
    if returncode:
        sys.exit(f"Error: 'samtools split' failed.\nStdout: {stdout}\nStderr: {stderr}\n")

    # convert readGroupId to filename friendly
    # only process the lanes output for given input bam
    for lane_bam in glob.glob(os.path.join(os.getcwd(), bam_base, "*.bam")):

        # remove file extension to get rg_id
        rg_id = os.path.splitext(os.path.basename(lane_bam))[0]

        # let's make sure RG_ID in lane bam exists in readgroup metadata, either matching read_group_id_in_bam or submitter_read_group_id
        # otherwise, it should be a lane that's expected to be ignored
        rg_id_found = False
        for rg in readgroups:
            if rg.get('file_r1') == os.path.basename(bam) and (rg.get('read_group_id_in_bam') == rg_id or
                    (not rg.get('read_group_id_in_bam') and rg['submitter_read_group_id'] == rg_id)):
                rg_id_found = True
                # rgs_with_lane_bam_produced.add(rg['submitter_read_group_id'])
                break

        if rg_id_found:
            rg_id_fn = readgroup_id_to_fname(rg_id, os.path.basename(bam), study_id, donor_id, sample_id)
            # 0x900 == 2304 == not primary alignment + supplementary alignments
            if rg['is_paired_end']:
              cmd = f"samtools collate -uO --threads {cpu} {lane_bam} | \
                samtools fastq -N -O -F 0x900 --threads {cpu} -0 {out_dir}/{rg_id_fn}_other.fq.gz \
              -1 {out_dir}/{rg_id_fn}_R1.fq.gz -2 {out_dir}/{rg_id_fn}_R2.fq.gz \
              -s {out_dir}/{rg_id_fn}_singleton.fq.gz - "
            else:
              cmd = f"samtools collate -uO --threads {cpu} {lane_bam} | \
              samtools fastq -N -O -F 0x900 --threads {cpu} -0 {out_dir}/{rg_id_fn}_other.fq.gz \
              -1 {out_dir}/{rg_id_fn}_R1.fq.gz -s {out_dir}/{rg_id_fn}_singleton.fq.gz - "
            
            stdout, stderr, returncode = run_cmd(cmd)
            if returncode:
              sys.exit(f"Error: 'samtools collate and fastq' failed.\nStdout: {stdout}\nStderr: {stderr}\n")
            
            sample_sheet[rg['submitter_read_group_id']] = {
              'file_r1': os.path.join(os.getcwd(), f'{out_dir}/{rg_id_fn}_R1.fq.gz'),
              'file_r2': os.path.join(os.getcwd(), f'{out_dir}/{rg_id_fn}_R2.fq.gz') if os.path.isfile(f'{out_dir}/{rg_id_fn}_R2.fq.gz') else 'No_File'
            }

            # retrieve read_group_info from metadata
            read_group_info = get_read_group_info(rg, study_id, donor_id, sample_id, specimen_id, specimen_type, tumour_normal_designation)
            
            if read_group_info:
              rg_kv = [ '@RG' ] + [ '%s:%s' % (k, v) for k, v in read_group_info.items() ]
              rg_array = "\'"+'\\t'.join(rg_kv)+"\'"

            sample_sheet[rg['submitter_read_group_id']].update({'read_group': rg_array}) 

        else:  # ignore lane bam without read group information in metadata, just produce a warning message here
            print("WARNING: Ignore lane BAM '%s' (split from input BAM '%s') that has no corresponding read group in the metadata" %
                  (lane_bam, os.path.basename(bam)), file=sys.stderr)
    for lane_bam in glob.glob(os.path.join(os.getcwd(), bam_base, "*.bam")):
      cmd = f"rm %s" % (lane_bam)
      stdout, stderr, returncode = run_cmd(cmd)
      if returncode:
        sys.exit(f"Error: 'Deleting file {lane_bam}' failed.\nStdout: {stdout}\nStderr: {stderr}\n")      

    return sample_sheet    

def bunzip2(fq_pair):
    bunzipped = []
    for f in fq_pair:
        if f and f.endswith('.bz2'):
            seq_file_bunzipped = os.path.join(
                                    os.environ.get("TMPDIR", ''),
                                    os.path.basename(re.sub(r'\.bz2$', '', f))
                                )
            cmd = 'bunzip2 -k -c %s > %s' % (f, seq_file_bunzipped)
            try:
                subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            except Exception as e:
                sys.exit("Unable to uncompress bz2 file: %s. Error: %s" % (f, e))

            bunzipped.append(seq_file_bunzipped)
        else:
            bunzipped.append(f)

    return bunzipped

def get_new_filename(fastq_old, rg_id_fn, r1_r2, outdir):
    if fastq_old.endswith('fq') or fastq_old.endswith('fastq'):
      fastq_new = os.path.join(os.getcwd(), outdir, f'{rg_id_fn}_{r1_r2}.fq')
    elif fastq_old.endswith('fq.gz') or fastq_old.endswith('fastq.gz'):
      fastq_new = os.path.join(os.getcwd(), outdir, f'{rg_id_fn}_{r1_r2}.fq.gz')
    else:
      sys.exit("Unsupported file format: %s." % fastq_old)
    
    return fastq_new
    
def get_read_group_info(read_group, study_id, donor_id, sample_id, specimen_id, specimen_type, tumour_normal_designation):

    read_group_info = {
        'ID': read_group['submitter_read_group_id'],
        'SM': sample_id,
        'LB': read_group['library_name'],
        'PU': read_group['platform_unit']
    }

    if read_group.get('insert_size'):
        read_group_info.update({'PI': read_group['insert_size']})
    if read_group.get('sample_barcode'):
        read_group_info.update({'BC': read_group['sample_barcode']})
    if read_group['experiment'].get('sequencing_center'):
        read_group_info.update({'CN': read_group['experiment']['sequencing_center']})
    if read_group['experiment'].get('platform'):
        read_group_info.update({'PL': read_group['experiment']['platform']})
    if read_group['experiment'].get('platform_model'):
        read_group_info.update({'PM': read_group['experiment']['platform_model']})
    if read_group['experiment'].get('sequencing_date'):
        read_group_info.update({'DT': read_group['experiment']['sequencing_date']})

    description = '|'.join([
                                read_group['experiment']['experimental_strategy'],
                                study_id,
                                specimen_id,
                                donor_id,
                                specimen_type,
                                tumour_normal_designation
                            ])

    read_group_info.update({'DS': description})

    return read_group_info

def main():
    """
    Python implementation of tool: prep-sample
    """

    parser = ArgumentParser(description='Tool: prep-sample')
    parser.add_argument("-s", "--input-files", dest="input_files", required=True,
                        help="Input files to process", type=str, nargs='+')
    parser.add_argument("-p", "--metadata-json", dest="metadata_json", required=True,
                        help="JSON file containing song analysis")
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, default="out",
                        help='Specify directory for output files')
    parser.add_argument("-n", "--cpus", type=int, default=cpu_count())

    args = parser.parse_args()

    if os.path.exists(args.outdir) and os.path.isdir(args.outdir):
      shutil.rmtree(args.outdir)
    os.mkdir(args.outdir)

    with open(args.metadata_json, 'r') as f:
      song_analysis = json.load(f)

    metadata_json = os.path.join(os.getcwd(), args.outdir, os.path.basename(args.metadata_json))
    os.symlink(os.path.abspath(args.metadata_json), metadata_json)
    study_id = song_analysis['studyId']
    donor_id = song_analysis['samples'][0]['donor']['donorId']
    sample_id = song_analysis['samples'][0]['sampleId']
    specimen_id = song_analysis['samples'][0]['specimen']['specimenId']
    if song_analysis['samples'][0]['donor']['gender'] == 'Female':
       sex = 'XX'  
    elif song_analysis['samples'][0]['donor']['gender'] == 'Male': 
       sex = 'XY'
    else:
       sex = 'NA'
    specimen_type = song_analysis['samples'][0]['specimen']['specimenType']
    tumour_normal_designation = song_analysis['samples'][0]['specimen']['tumourNormalDesignation']
    status = '0' if tumour_normal_designation == 'Normal' else '1'
    
    if song_analysis.get('workflow'):
      genome_build = song_analysis['workflow']["genome_build"]
    else:
      genome_build = None

    analysis_type = song_analysis['analysisType']['name']
    output_sample_sheet = f'{args.outdir}/{sample_id}_{analysis_type}_sample_sheet.csv'
    experiment=song_analysis['experiment']['experimental_strategy']


    sample_sheet = dict()
    if analysis_type == 'sequencing_experiment':
      read_group_count = song_analysis['read_group_count']
      filepair_map_to_readgroup = group_readgroup_by_filepair(song_analysis)

      for fp in filepair_map_to_readgroup:
        if filepair_map_to_readgroup[fp]['format'] == 'BAM':
          # for bam just need fp[0] since fp[1] is either the same as fp[0] or None
          sample_sheet = generate_fastqs_from_bam(filename_to_file(fp, args.input_files)[0],
                                      filepair_map_to_readgroup[fp]['read_groups'],
                                      args.cpus, sample_sheet, study_id, donor_id, sample_id, 
                                      args.outdir, specimen_id, specimen_type, 
                                      tumour_normal_designation)
        else: # FASTQ must be one read group
          fq_pair = filename_to_file(fp, args.input_files)
          fastq_pair = bunzip2(fq_pair)
          rg = filepair_map_to_readgroup[fp]['read_groups'][0]
          rg_id = rg['submitter_read_group_id']
          rg_id_fn = readgroup_id_to_fname(rg_id, '', study_id, donor_id, sample_id)
          file_r1_new = get_new_filename(fastq_pair[0], rg_id_fn, "R1", args.outdir)
          os.symlink(os.path.abspath(fastq_pair[0]), file_r1_new)
          if fastq_pair[1]:
            file_r2_new = get_new_filename(fastq_pair[1], rg_id_fn, "R2", args.outdir)
            os.symlink(os.path.abspath(fastq_pair[1]), file_r2_new)
          else:
            file_r2_new = 'No_File'

          sample_sheet[rg_id] = {
            'file_r1': file_r1_new,
            'file_r2': file_r2_new
          }
        
          # retrieve read_group_info from metadata
          read_group_info = get_read_group_info(rg, study_id, donor_id, sample_id, specimen_id, specimen_type, tumour_normal_designation)

          if read_group_info:
            rg_kv = [ '@RG' ] + [ '%s:%s' % (k, v) for k, v in read_group_info.items() ]
            rg_array = "\'"+'\\t'.join(rg_kv)+"\'"

          sample_sheet[rg_id].update({'read_group': rg_array}) 
        
          
      # now we check whether all read groups in metadata have produced lane fastq
      rgs_missed_lane = set()
      for rg in song_analysis['read_groups']:
          if rg['submitter_read_group_id'] not in sample_sheet:
              rgs_missed_lane.add(rg['submitter_read_group_id'])

      if rgs_missed_lane:  # throw error here if that happens
          sys.exit("Error: no lane BAM has been generated for some read groups: '%s'. "
                  "Please make sure supplied sequencing files and metadata are correct." % "', '".join(rgs_missed_lane))

      with open(output_sample_sheet, 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerow(['analysis_type','study_id','patient','sex','status','sample','lane','fastq_1','fastq_2','read_group','single_end','read_group_count',"experiment", 'analysis_json'])
        for k,v in sample_sheet.items():
          single_end = True if v['file_r2'] == 'No_File' else False
          csvwriter.writerow([analysis_type, study_id, donor_id, sex, status, sample_id, k, v['file_r1'], v['file_r2'], v['read_group'], single_end, read_group_count,experiment, metadata_json])
    
    elif analysis_type == 'sequencing_alignment':
      for fp in args.input_files:
        if fp.endswith('cram'): 
          cram = os.path.join(os.getcwd(), args.outdir, os.path.basename(fp))
          os.symlink(os.path.abspath(fp), cram)
        elif fp.endswith('crai'):
          crai = os.path.join(os.getcwd(), args.outdir, os.path.basename(fp))
          os.symlink(os.path.abspath(fp), crai)
        else:
          sys.exit("Error: not supported input file format")
      with open(output_sample_sheet, 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerow(['analysis_type','study_id','patient','sex','status','sample','cram','crai',"genome_build",'experiment', 'analysis_json'])
        csvwriter.writerow([analysis_type, study_id, donor_id, sex, status, sample_id, cram, crai, genome_build,experiment, metadata_json])

    elif analysis_type == 'variant_calling':
      for fp in song_analysis['files']:
         if not fp['fileType'] == 'VCF': continue
         variantcaller = fp['info']['analysis_tools'][0]
      for fp in args.input_files:
        if fp.endswith('vcf.gz'): 
          vcf = os.path.join(os.getcwd(), args.outdir, os.path.basename(fp))
          os.symlink(os.path.abspath(fp), vcf)
        elif fp.endswith('vcf.gz.tbi'):
          tbi = os.path.join(os.getcwd(), args.outdir, os.path.basename(fp))
          os.symlink(os.path.abspath(fp), tbi)
        else:
          sys.exit("Error: not supported input file format")
      with open(output_sample_sheet, 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerow(['analysis_type','study_id','patient','sex','sample','variantcaller','vcf','tbi',"genome_build",'experiment', 'analysis_json'])
        csvwriter.writerow([analysis_type, study_id, donor_id, sex, sample_id, variantcaller, vcf, tbi ,genome_build,experiment, metadata_json])  

    elif analysis_type == 'qc_metrics':
      with open(output_sample_sheet, 'w', newline='') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerow(['analysis_type','study_id','patient','sex','status','sample','qc_tools','qc_file',"genome_build", 'experiment','analysis_json'])

        for fp in args.input_files:
          for fq in song_analysis['files']:
            if not fq.get('fileName') == os.path.basename(fp): continue
            qc_file = os.path.join(os.getcwd(), args.outdir, os.path.basename(fp))
            os.symlink(os.path.abspath(fp), qc_file)
            qc_tools = ','.join(fq['info']['analysis_tools'])

          csvwriter.writerow([analysis_type, study_id, donor_id, sex, status, sample_id, qc_tools, qc_file, genome_build, experiment, metadata_json]) 

if __name__ == "__main__":
    main()

