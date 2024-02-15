#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".bam",
        ".cram",
    )

    def __init__(
        self,
        #sample_col="sample",
        #first_col="bam_cram",
        analysis_type_col = 'analysis_type',
        study_id_col = 'study_id',
        patient_col = 'patient',
        sex_col = 'sex',
        status_col = 'status',
        sample_col = 'sample',
        cram_col = 'cram',
        crai_col = 'crai',
        experiment_col = 'experiment',
        genome_build_col = "genome_build",
        analysis_json_col = 'analysis_json',
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.
analysis_type,study_id,patient,sex,status,sample,lane,fastq_1,fastq_2,read_group,single_end,read_group_count,analysis_json
        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._analysis_type_col = analysis_type_col
        self._study_id_col = study_id_col
        self._patient_col = patient_col
        self._sex_col = sex_col
        self._status_col = status_col
        self._sample_col = sample_col
        self._cram_col = cram_col
        self._crai_col = crai_col
        self._experiment_col = experiment_col
        self._genome_build_col = genome_build_col
        self._analysis_json_col = analysis_json_col
        self._seen = []
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        #{"analysis_type","study_id","patient","sex","status","sample","cram","crai","analysis_json"}
        self._validate_analysis_type(row) if row.get(self._analysis_type_col) else ""
        self._validate_sex(row) if row.get(self._sex_col) else ""
        self._validate_study_id(row) if row.get(self._study_id_col) else ""
        self._validate_patient(row) if row.get(self._patient_col) else ""
        self._validate_status(row) if row.get(self._status_col) else ""
        self._validate_sample(row)
        self._validate_cram(row)
        self._validate_crai(row)
        self._validate_experiment(row) if row.get(self._experiment_col) else ""
        self._validate_genome_build(row) if row.get(self._genome_build_col) else ""
        self._validate_analysis_json(row) if row.get(self._analysis_json_col) else ""


        tmp_dict={
            "analysis_type" : row[self._analysis_type_col] if row.get(self._analysis_type_col) else "sequencing_alignment",
            "study_id" : row[self._study_id_col] if row.get(self._study_id_col) else "LOCAL",
            "patient" : row[self._patient_col] if row.get(self._patient_col) else row[self._sample_col],
            "sex" : row[self._sex_col] if row.get(self._sex_col) else "NA",
            "status" : row[self._status_col] if row.get(self._status_col) else "0",
            "sample" : row[self._sample_col],
            "cram" : row[self._cram_col],
            "crai" : row[self._crai_col],
            "experiment": row[self._experiment_col] if row.get(self._experiment_col) else "WGS",
            "genome_build": row[self._genome_build_col] if row.get(self._genome_build_col) else "GRCh38",
            "analysis_json": row[self._analysis_json_col] if row.get(self._analysis_json_col) else None
            }

        self._seen.append(row)
        self.modified.append(tmp_dict)


    def _validate_analysis_type(self, row):
        """Assert that expected analysis is correct."""
        if len(row[self._analysis_type_col]) <= 0:
            raise AssertionError("'analysis_type' input is required.")
        if row[self._analysis_type_col]!="sequencing_alignment":
            raise AssertionError("analysis_type for \"DNA Alignment QC\" should be  \"sequencing_alignment\"")

    def _validate_study_id(self, row):
        """Assert that expected study_id is correct."""
        if len(row[self._study_id_col]) <= 0:
            raise AssertionError("'study_id' input is required.")

    def _validate_patient(self, row):
        """Assert that expected patient is correct."""
        if len(row[self._patient_col]) <= 0:
            raise AssertionError("'patient' input is required.")

    def _validate_sex(self, row):
        """Assert that expected sex is correct."""
        if len(row[self._sex_col]) <= 0:
            raise AssertionError("'analysis_type' input is required.")
        if row[self._sex_col]!="XX" and row[self._sex_col]!="XY" and row[self._sex_col]!="NA":
            raise AssertionError("sex should be one of the following values : XX,XY,NA")

    def _validate_status(self, row):
        """Assert that expected tumour status is correct."""
        if len(row[self._status_col]) <= 0:
            raise AssertionError("'status' input is required.")
        if row[self._status_col]!="1" and row[self._status_col]!="0":
            raise AssertionError("Tumour status should be \"0\" is normal else \"1\"")

    def _validate_sample(self, row):
        """Assert that expected sample is correct."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("'sample' input is required.")
    
    def _validate_cram(self, row):
        """Assert that expected cram is correct."""
        if len(row[self._cram_col]) <= 0:
            raise AssertionError("'cram' input is required.")
        if not row[self._cram_col].endswith(".cram"):
            raise AssertionError("'cram' input format is incorrect, ensure file ends with '.cram'")

    def _validate_crai(self, row):
        """Assert that expected crai is correct."""
        if len(row[self._crai_col]) <= 0:
            raise AssertionError("'crai' input is required.")
        if not row[self._crai_col].endswith(".crai"):
            raise AssertionError("'crai' input format is incorrect, ensure file ends with '.crai'")
        if row[self._crai_col].split("/")[-1].replace(".cram.crai","")!=row[self._cram_col].split("/")[-1].replace(".cram",""):
            raise AssertionError("'cram' and 'crai' file name bodies do not match.")

    def _validate_experiment(self, row):
        """Assert that expected Experiment is correct."""
        if len(row[self._experiment_col]) <= 0:
            raise AssertionError("'experiment' input is required.")
        for val in ["WGS","WXS","RNA-Seq","Bisulfite-Seq","ChIP-Seq","Targeted-Seq"]:
            if val==row[self._experiment_col]:
                return
        raise AssertionError("'experiment' type does not match the following: \"WGS\",\"WXS\",\"RNA-Seq\",\"Bisulfite-Seq\",\"ChIP-Seq\",\"Targeted-Seq\".")


    def _validate_analysis_json(self, row):
        """Assert that expected analysis_json is correct."""
        if len(row[self._analysis_json_col]) <= 0:
            raise AssertionError("'analysis_json' input is required.")
        if not row[self._analysis_json_col].endswith(".json"):
            raise AssertionError("'analysis_json' input should have the suffix \".json\".")

    def _validate_genome_build(self, row):
        """Assert that expected genome_build is correct."""
        if len(row[self._genome_build_col]) <= 0:
            raise AssertionError("'genome_build' input is required.")

    def validate_unique_values(self,col):
        """
        Assert a single unique value exists in array
        """
        if len(set([z[col] for z in self.modified]))!=len([z[col] for z in self.modified]):
                raise AssertionError("Errors duplicates values detected for '%s'. Each row should have an unique value" % (col))
                sys.exit(1)

    def validate_common_values(self,col):
        """
        Assert each value in array is unique
        """
        if len(set([z[col] for z in self.modified]))!=1:
            raise AssertionError("Errors multiple values detected for '%s'. Only a single value should be used" % (col))
            sys.exit(1)



def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
    file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
    file_out (pathlib.Path): Where the validated and transformed samplesheet should
        be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,

    analysis_type,study_id,patient,sex,status,sample,cram,crai,genome_build,analysis_json
    sequencing_alignment,TEST-QA,DO262466,XY,1,SA622744,TEST-QA.DO262466.SA622744.wxs.20210712.aln.cram,TEST-QA.DO262466.SA622744.wxs.20210712.aln.cram.crai,WXS,hg38,4f6d6ddf-3759-4a30-ad6d-df37591a3033.analysis.json
    """
    required_columns = {"sample","cram","crai"}
    conditional_columns = {"study_id","sex","patient","status","experiment","analysis_json"}

    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames) and not conditional_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)

        for col in["sample","study_id","sex","patient","experiment","status","analysis_json"]:
            checker.validate_common_values(col)
        for col in ["cram","crai"]:
            checker.validate_unique_values(col)

    header = checker.modified[0].keys()
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog=\
        '''
Check that the tabular samplesheet has the structure expected by nf-core pipelines.

Validate the general shape of the table, expected columns, and each row. Also add
an additional column which records whether one or two FASTQ reads were found.

Args:
file_in (pathlib.Path): The given tabular samplesheet. The format can be either
        CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
file_out (pathlib.Path): Where the validated and transformed samplesheet should
    be created; always in CSV format.

Example:
    This function checks that the samplesheet follows the following structure,

    analysis_type,study_id,patient,sex,status,sample,cram,crai,analysis_json
    sequencing_alignment,TEST-QA,DO262466,XY,1,SA622744,TEST-QA.DO262466.SA622744.wxs.20210712.aln.cram,TEST-QA.DO262466.SA622744.wxs.20210712.aln.cram.crai,4f6d6ddf-3759-4a30-ad6d-df37591a3033.analysis.json
    ''',

        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())