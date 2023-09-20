## Introduction

**icgc-argo-workflows/dnaalnqc** is a reproducible bioinformatics workflow that can be used to obtain QC metrics from tumour/normal paired, tumour-only, normal-only WGS/WXS/Targeted-Seq aligned reads. It has been created to support quality control efforts within ICGC-ARGO project. The aggregated QC metrics are formed to align with the [GA4GH WGS_Quality_Control_Standards](https://www.ga4gh.org/product/wgs-quality-control-standards/).  

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

The workflow is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. 
The workflow has adopted [nf-core](https://nf-co.re/) framework and best practice guidelines to ensure reproducibility, portability and scalability. Where possible, many processes have been installed from [nf-core/modules](https://github.com/nf-core/modules). Moreover, ICGC ARGO specific modules have been installed form [icgc-argo-workflows/argo-modules](https://github.com/icgc-argo-workflows/argo-modules), which hosts ARGO reusable modules across all ICGC ARGO pipelines!

## Requirements

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Stage the required [reference files](#references) 

## Quick start
1. Test the workflow running in `Local` mode on a minimal dataset with a single command:

   ```bash
   nextflow run icgc-argo-workflows/dnaalnqc \
     -profile test,standard \
     --outdir <OUTDIR>
   ```

2. Test the workflow running in `RDPC` mode with a single command if you have access to `RDPC-QA` env and have your valid api_token available:

   ```bash
   nextflow run icgc-argo-workflows/dnaalnqc \
     -profile test_rdpc_qa,standard \
     --api_token <YOUR_API_TOKEN> \
     --reference_base <REFERENCE_BASE> \
     --outdir <OUTDIR>
   ```

## Usage

### Workflow summary
Depending on where the input data are coming from and output data are sending to, the workflow can be running in two modes: `Local` and `RDPC` . The major tasks performed in the workflow are:
- (`RDPC` mode only) Download input sequencing metadata/data from data center using SONG/SCORE client tools
- Perform `Samtools Stats` to collect metrics related to sequencing library quality and read characteristics
- Perform `Picard CollectOxoGMetrics` to collect OxoG metrics
- Perform `Picard CollectQualityYieldMetrics` to collect read characteristics metrics
- Perform `Picard CollectWGSMetrics` to collect coverage metrics for WGS
- Perform `Picard CollectHSMetrics` to collect coverage metrics for WXS and Targeted-Seq
- Perfomr `GATK CalculateContamination` to collect cross sample contamination rate
- Perform `MultiQC` analysis to generate aggregated results
- (`RDPC` mode only) Generate SONG metadata for all collected QC metrics files and upload QC files to SONG/SCORE


### References
- Reference genome: 
  - GRCh38 reference genome fasta file. The file can be downloaded by:
    ```bash
    wget https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.fa
    ``` 

  - GRCh38 reference genome fasta index file. The file can be downloaded by:
    ```bash
    wget https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.fa.fai
    ```

  - GRCh38 reference genome sequence dictionary file. The file can be downloaded by:
    ```bash
    wget https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.dict
    ```
- GATK resources: 
  - `germline_resource` and index files. The files can be downloaded by:
    ```bash
    wget https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/gatk-resources/af-only-gnomad.pass-only.biallelic.snp.hg38.vcf.gz
    wget https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/gatk-resources/af-only-gnomad.pass-only.biallelic.snp.hg38.vcf.gz.tbi
    ``` 
- Autosome non-gap regions
  - `autosome_non_gap` bed file was downloaded from [NPM-sample-qc](https://raw.githubusercontent.com/c-BIG/NPM-sample-qc/master/resources/autosomes_non_gap_regions.bed) and staged under project folder [assets](https://github.com/icgc-argo-workflows/dnaalnqc/tree/main/assets)


### Inputs
#### Local mode
First, prepare a samplesheet with your input data that looks as following:

`samplesheet.csv`:

```csv
biosample_id,bam_cram,patient(optional),status(optional),sex(optional)
CONTROL_REP1_SAMPLE0,CONTROL_REP_0.bam,CONTROL_REP1_DONOR,0,XX
CONTROL_REP1_SAMPLE1,CONTROL_REP_1.bam,CONTROL_REP1_DONOR,1,XX
```

Each row represents an aligned BAM or CRAM from a sample.

Then, you need to download all required [reference files](#references), and stage them into a directory <REFERENCE_BASE>

Now, you can run the workflow using:

```bash
nextflow run icgc-argo-workflows/dnaalnqc \
   -profile resource,<docker/singularity> \
   --local_mode true \
   --input samplesheet.csv \
   --reference_base <REFERENCE_BASE> \
   --outdir <OUTDIR>
```

#### RDPC mode
You can run the workflow in RDPC mode by using:
```bash
nextflow run icgc-argo-workflows/dnaalnqc \
  -profile <rdpc,rdpc_qa,rdpc_dev>,<docker/singularity> \
  --local_mode false \
  --study_id <STUDY_ID> \
  --analysis_ids <ANALYSIS_IDS> \
  --api_token <YOUR_API_TOKEN> \ 
  --reference_base <REFERENCE_BASE> \
  --outdir <OUTDIR>
```

> **NOTE**
> Please provide workflow parameters via the CLI or Nextflow `-params-file` option. 

### Outputs
Upon completion, you can find the aggregated QC metrics under directory:
```
/path/to/outdir/metrics/<sample_id>.multiqc_data.json
```

## Credits

icgc-argo-workflows/dnaalnqc was mostly written by Linda Xiang (@lindaxiang), with contributions from 
Andrej Benjak, Charlotte Ng, Desiree Schnidrig, Edmund Su, Miguel Vazquez, Morgan Taschuk, Raquel Manzano Garcia, Romina Royo and ICGC-ARGO Quality Control Working Group.  

Authors (alphabetical)
- Andrej Benjak
- Charlotte Ng
- Desiree Schnidrig
- Edmund Su
- Linda Xiang
- Miguel Vazquez
- Morgan Taschuk
- Raquel Manzano Garcia
- Romina Royo

## Citations
<!-- If you use  icgc-argo-workflows/dnaalnqc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
