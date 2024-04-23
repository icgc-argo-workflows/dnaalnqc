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
        analysis_type_col = 'analysis_type',
        study_id_col = 'study_id',
        patient_col = 'patient',
        sex_col = 'sex',
        status_col = 'status',
        sample_col = 'sample',
        lane_col = 'lane',
        fastq_1_col = 'fastq_1',
        fastq_2_col = 'fastq_2',
        library_name_col = 'library_name',
        platform_unit_col = 'platform_unit',
        platform_col = 'platform',
        sequencing_center_col = 'sequencing_center',
        sequencing_date_col = 'sequencing_date',
        platform_model_col = 'platform_model',
        single_end_col = 'single_end',
        read_group_count_col = 'read_group_count',
        experiment_col = 'experiment',
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
        self._lane_col = lane_col
        self._fastq_1_col = fastq_1_col
        self._fastq_2_col = fastq_2_col
        self._library_name_col = library_name_col
        self._platform_unit_col = platform_unit_col
        self._platform_col = platform_col
        self._sequencing_center_col = sequencing_center_col
        self._sequencing_date_col = sequencing_date_col
        self._platform_model_col = platform_model_col
        self._single_end_col = single_end_col
        self._read_group_count_col = read_group_count_col
        self._experiment_col = experiment_col
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
        self._validate_analysis_type(row) if row.get(self._analysis_type_col) else ""
        self._validate_sex(row) if row.get(self._sex_col) else ""
        self._validate_study_id(row) if row.get(self._study_id_col) else ""
        self._validate_patient(row) if row.get(self._patient_col) else ""
        self._validate_sex(row) if row.get(self._sex_col) else ""
        self._validate_status(row) if row.get(self._status_col) else ""
        self._validate_sample(row)
        self._validate_lane(row)
        self._validate_single_end(row)
        self._validate_fastq_1(row)
        self._validate_fastq_2(row)
        self._validate_library_name(row)
        self._validate_platform_unit(row)
        self._validate_platform_col(row) if row.get(self._platform_col) else ""
        self._validate_sequencing_center_col(row) if row.get(self._sequencing_center_col) else ""
        self._validate_sequencing_date_col(row) if row.get(self._sequencing_date_col) else ""
        self._validate_platform_model_col(row) if row.get(self._platform_model_col) else ""
        self._validate_read_group_count(row)
        self._validate_experiment(row) if row.get(self._experiment_col) else ""
        self._validate_analysis_json(row) if row.get(self._analysis_json_col) else ""

        tmp_dict={
            "analysis_type" : row[self._analysis_type_col] if row.get(self._analysis_type_col) else "sequencing_experiment",
            "study_id" : row[self._study_id_col] if row.get(self._study_id_col) else "LOCAL",
            "patient" : row[self._patient_col] if row.get(self._patient_col) else row[self._sample_col],
            "sex" : row[self._sex_col] if row.get(self._sex_col) else "NA",
            "status" : row[self._status_col] if row.get(self._status_col) else "0",
            "sample" : row[self._sample_col],
            "lane" : row[self._lane_col],
            "fastq_1" : row[self._fastq_1_col],
            "fastq_2" : row[self._fastq_2_col] if row.get(self._fastq_2_col) else "NO_FILE",
            "single_end" : row[self._single_end_col].lower(),
            "read_group_count" : row[self._read_group_count_col],
            "experiment" : row[self._experiment_col] if row.get(self._experiment_col) else "WGS",
            "analysis_json": row[self._analysis_json_col] if row.get(self._analysis_json_col) else None
            }

        read_group_info=[]
        description=[]

        for col in [
            'experiment',
            'study_id',
            'experiment',
            'patient',
            'sample',
            'status'
        ]:
            if tmp_dict.get(col):
                if col=='status':
                    if tmp_dict['status']==1:
                        description.append("Tumour")
                    else:
                        description.append("Normal")
                        continue
                description.append(tmp_dict[col])

        for col,id in zip(
            [
                self._lane_col,
                self._sample_col,
                self._library_name_col,
                self._platform_unit_col,
                self._sequencing_center_col,
                self._platform_col,
                self._platform_model_col,
                self._sequencing_date_col
            ],
            ["ID","SM","LB","PU","CN","PL","PM","DT"]):
            if row.get(col):
                read_group_info.append("%s:%s" % (id,row[col]))

        tmp_dict['read_group']="'@RG\\t%s\\tDS:%s'" % ("\\t".join(read_group_info),"|".join(description))

        self._seen.append(row)
        self.modified.append(tmp_dict)

    def _validate_analysis_type(self, row):
        """Assert that expected analysis is correct."""
        if len(row[self._analysis_type_col]) <= 0:
            raise AssertionError("'analysis_type' input is required.")
        if row[self._analysis_type_col]!="sequencing_experiment":
            raise AssertionError("analysis_type for \"DNA Alignment\" should be  \"sequencing_experiment\"")

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
    

    def _validate_lane(self, row):
        """Assert that expected lane is correct."""
        if len(row[self._lane_col]) <= 0:
            raise AssertionError("'lane' input is required.")
    

    def _validate_fastq_1(self, row):
        """Assert that expected fastq_1 is correct."""
        if len(row[self._fastq_1_col]) <= 0:
            raise AssertionError("'fastq_1' input is required.")
        if not (
            row[self._fastq_1_col].endswith(".fq.gz") or 
            row[self._fastq_1_col].endswith(".fastq.gz") or
            row[self._fastq_1_col].endswith(".bam")
            ):
            raise AssertionError("'fastq_1' incorrect format detected.")
    

    def _validate_fastq_2(self, row):
        """Assert that expected fastq_2 is correct."""
        if row[self._single_end_col].lower()=="true":
            return 

        if len(row[self._fastq_2_col]) <= 0:
            raise AssertionError("'fastq_2' input is required.")
        if row[self._fastq_2_col].endswith(".fastq.gz"):
            if row[self._fastq_2_col].split("/")[-1].replace("R2.fastq.gz","").replace("r2.fastq.gz","")!=row[self._fastq_1_col].split("/")[-1].replace("R1.fastq.gz","").replace("r1.fastq.gz",""):
                raise AssertionError("'fastq_1' and 'fastq_2' prefix differ.")
        if row[self._fastq_2_col].endswith(".fq.gz"):
            if row[self._fastq_2_col].split("/")[-1].replace("R2.fq.gz","").replace("r2.fq.gz","")!=row[self._fastq_1_col].split("/")[-1].replace("R1.fq.gz","").replace("r1.fq.gz",""):
                raise AssertionError("'fastq_1' and 'fastq_2' prefix differ.")
        if row[self._fastq_2_col].endswith(".bam"):
            if row[self._fastq_2_col]!=row[self._fastq_1_col]:
                raise AssertionError("'fastq_1' and 'fastq_2' prefix differ.")

    def _validate_single_end(self, row):
        """Assert that expected single_end is correct."""
        if len(row[self._single_end_col]) <= 0:
            raise AssertionError("'single_end' input is required.")
        if row[self._single_end_col].lower()!="true" and row[self._single_end_col].lower()!="false":
            raise AssertionError("'single_end' should be specifed as \"True\" or \"False\".") 
    

    def _validate_read_group_count(self, row):
        """Assert that expected read_group_count is correct."""
        if len(row[self._read_group_count_col]) <= 0:
            raise AssertionError("'read_group_count' input is required.")
    
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

    def _validate_library_name(self, row):
        """Assert that expected library_name is correct."""
        if len(row[self._library_name_col]) <= 0:
            raise AssertionError("'library_name' input is required.")

    def _validate_platform_unit(self, row):
        """Assert that expected platform_unit is correct."""
        if len(row[self._platform_unit_col]) <= 0:
            raise AssertionError("'platform_unit' input is required.")

    def _validate_platform_col(self, row):
        """Assert that expected platform is correct."""
        if len(row[self._platform_col]) <= 0:
            raise AssertionError("'platform' input is required.")

    def _validate_sequencing_center_col(self, row):
        """Assert that expected sequencing_center is correct."""
        if len(row[self._sequencing_center_col]) <= 0:
            raise AssertionError("'sequencing_center' input is required.")

    def _validate_sequencing_date_col(self, row):
        """Assert that expected sequencing_date is correct."""
        if len(row[self._sequencing_date_col]) <= 0:
            raise AssertionError("'sequencing_date' input is required.")

    def _validate_platform_model_col(self, row):
        """Assert that expected platform_model is correct."""
        if len(row[self._platform_model_col]) <= 0:
            raise AssertionError("'platform_model' input is required.")

    def validate_unique_fastq(self):
        """
        Assert that the combination of FASTQ filename is unique.
        """
        tmp=[z['fastq_1'] for z in self.modified]+[z['fastq_2'] for z in self.modified]

        for iter in range(0,len(tmp)):
            current_val=tmp.pop(0)
            if current_val.endswith(".fastq.gz"):
                continue
            if current_val.endswith(".fq.gz"):
                continue
            if current_val=='NO_FILE':
                continue
            if current_val in tmp:
                raise AssertionError("Errors multiple instances of file '%s' detected" % (current_val))
                sys.exit(1)
            else:
                raise AssertionError("Unexpected file format detected for '%s'" % (current_val))

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
    required_columns = {"sample","lane","fastq_1","fastq_2","single_end","read_group_count","library_name","platform_unit"}
    conditional_columns = {"study_id","sex","patient","status","experiment","analysis_json","platform","sequencing_center","sequencing_date","platform_model"}

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
        checker.validate_unique_fastq()
        for col in["sample","study_id","sex","patient","experiment","read_group_count","status","analysis_json"]:
            checker.validate_common_values(col)
        for col in ["lane"]:
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

    analysis_type,study_id,patient,sex,status,sample,lane,fastq_1,fastq_2,read_group,single_end,read_group_count,analysis_json
    sequencing_experiment,TEST-QA,DO263089,XX,1,SA624380,C0HVY.2,TEST-QA.DO263089.SA624380.C0HVY.2.8775eee1cacedc27428856591023d837_R1.fq.gz,TEST-QA.DO263089.SA624380.C0HVY.2.8775eee1cacedc27428856591023d837_R2.fq.gz,'@RG\\tID:C0HVY.2\\tSM:SA624380\\tLB:Pond-147580\\tPU:74_8a\\tPI:298\\tCN:EXT\\tPL:ILLUMINA\\tPM:HiSeq 2000\\tDT:2014-12-12\\tDS:WGS|TEST-QA|SP224367|DO263089|Cell line - derived from tumour|Tumour',False,3,WXS,875ef550-e536-4456-9ef5-50e5362456df.analysis.json
    sequencing_experiment,TEST-QA,DO263089,XX,1,SA624380,D0RE2.1,TEST-QA.DO263089.SA624380.D0RE2.1.b8ac1a3b5b52ced6068b28c4e9b4e5e9_R1.fq.gz,TEST-QA.DO263089.SA624380.D0RE2.1.b8ac1a3b5b52ced6068b28c4e9b4e5e9_R2.fq.gz,'@RG\\tID:D0RE2.1\\tSM:SA624380\\tLB:Pond-147580\\tPU:74_8b\\tPI:298\\tCN:EXT\\tPL:ILLUMINA\\tPM:HiSeq 2000\\tDT:2014-12-12\\tDS:WGS|TEST-QA|SP224367|DO263089|Cell line - derived from tumour|Tumour',False,3,WXS,875ef550-e536-4456-9ef5-50e5362456df.analysis.json
    sequencing_experiment,TEST-QA,DO263089,XX,1,SA624380,D0RH0.2,TEST-QA.DO263089.SA624380.D0RH0.2.231146e66d802729c719428e33e555a8_R1.fq.gz,TEST-QA.DO263089.SA624380.D0RH0.2.231146e66d802729c719428e33e555a8_R2.fq.gz,'@RG\\tID:D0RH0.2\\tSM:SA624380\\tLB:Pond-147580\\tPU:74_8c\\tPI:298\\tCN:EXT\\tPL:ILLUMINA\\tPM:HiSeq 2000\\tDT:2014-12-12\\tDS:WGS|TEST-QA|SP224367|DO263089|Cell line - derived from tumour|Tumour',False,3,WXS,875ef550-e536-4456-9ef5-50e5362456df.analysis.json
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