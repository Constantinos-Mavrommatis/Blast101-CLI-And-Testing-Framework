import argparse
import os.path
import re
import programme_settings
from Bio import SeqIO
import sys


programme_settings.read()

#Default settings
default_database = programme_settings.settings["DEFAULT"]["database"]
default_query = programme_settings.settings["DEFAULT"]["query_sequence"]
default_blosum = (programme_settings.settings["DEFAULT"]["blosum"])
default_gap = (programme_settings.settings["DEFAULT"]["seq_gap"])

#Blast settings
default_blast_max_alignments = (programme_settings.settings["BLAST"]["max_alignments"])
default_blast_max_score = (programme_settings.settings["BLAST"]["max_scores"])
default_blast_word_size = (programme_settings.settings["BLAST"]["word_size"])

#SW settings
default_sw_max_score = (programme_settings.settings["SWSEARCH"]["max_sw_scores"])

def validate_sw_max_score(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer.")

    max_allowed = 5000
    if ivalue <= 0 or ivalue > max_allowed:
        raise argparse.ArgumentTypeError(f"Value must be less than or equal to {max_allowed} but bigger than 0.")

    return ivalue

def validate_word_size(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer.")

    min_allowed, max_allowed = 2, 6
    if not (min_allowed <= ivalue <= max_allowed):
        raise argparse.ArgumentTypeError(f"Must be between {min_allowed} and {max_allowed}.")

    return ivalue

def validate_blast_maximum_score(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer.")

    min_allowed, max_allowed = 1, 5000
    if not (min_allowed <= ivalue <= max_allowed):
        raise argparse.ArgumentTypeError(f"Gap penalty must be between {min_allowed} and {max_allowed}.")

    return ivalue

def validate_gap(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer.")

    min_allowed, max_allowed = -20, -1
    if not (min_allowed <= ivalue <= max_allowed):
        raise argparse.ArgumentTypeError(f"Gap penalty must be between {min_allowed} and {max_allowed}.")

    return ivalue

def validate_max_alignments(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer.")

    max_allowed = 1000
    if ivalue > max_allowed:
        raise argparse.ArgumentTypeError(f"Value must be less than or equal to {max_allowed}.")

    return ivalue

def validate_blosum(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer.")

    allowed_values = [45, 50, 62, 80, 90]
    if ivalue not in allowed_values:
        raise argparse.ArgumentTypeError(
            f"Must be one of {allowed_values}, but got {ivalue}."
        )
    return ivalue

def check_query_file(query_path):
    """
    """
    # 1. Check file existence
    if not os.path.isfile(query_path):
        print(f"Error: Query file '{query_path}' does not exist!")
        sys.exit(1)

    # 2. Check file size
    if os.stat(query_path).st_size == 0:
        print(f"Error: Query file '{query_path}' is empty!")
        sys.exit(1)

    # 4. Parse the first record with Biopython
    with open(query_path, "r") as fh:
        record = next(SeqIO.parse(fh, "fasta"))

    if re.fullmatch(r"[ATGC]+", str(record.seq.upper())):
        print(f"Error: Query file {query_path} is not a protein sequence!")
        sys.exit(1)

    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBZXOJU]+", str(record.seq.upper())):
        wrong_character = re.findall("[^ACDEFGHIKLMNPQRSTVWYBZXUOJ]", str(record.seq.upper()))
        print(f"Error: Record '{record.id}' in '{query_path}' has a wrong character '{wrong_character}'")
        sys.exit(1)

    return record.seq

def check_db_file(db_path):
    """
    Validates that 'db_path' is:
      1) An existing, non-empty file.
      2) A valid FASTA where *each* record starts with '>' and
         is strictly A/T/G/C if you need that restriction.
    Exits the program if any check fails.
    """
    # 1. Check file existence
    if not os.path.isfile(db_path):
        print(f"Error: Database file '{db_path}' does not exist!")
        sys.exit(1)

    # 2. Check file size
    if os.stat(db_path).st_size == 0:
        print(f"Error: Database file '{db_path}' is empty!")
        sys.exit(1)

    # 4. Parse with Biopython, check each record
    with open(db_path, "r") as seq:
        for record in SeqIO.parse(seq, "fasta"):
            seq_upper = record.seq.upper()

        #if re.fullmatch(r"[ATGC]+", str(seq_upper)):
            #print(f"Error: Database file {db_path} has a DNA sequence at '{record.id}'")
            #sys.exit(1)

            if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBZXOJU]+", str(seq_upper)):
                wrong_character=re.findall("[^ACDEFGHIKLMNPQRSTVWYBZXUOJ]", str(seq_upper))
                print(f"Error: Record '{record.id}' in '{db_path}' has a wrong character '{wrong_character}'")
                sys.exit(1)
    return db_path

def parse_args():

    parser = argparse.ArgumentParser(prog="Blast 101 Search BETA version Copyright (c) 2025 Simon Tomlinson University of Edinburgh",
                                     usage="[--db <file>] [--query <file>] [--method blast|sw|stats|tests]",
                                     description="Blast 101 Search: A tool for running BLAST and Smith-Waterman "
                                                 "searches, with default values for arguments read from "
                                                 "'settings.ini'. If no command-line arguments are specified, "
                                                 "the default BLAST run uses the values found in 'settings.ini'.",
                                     epilog="LICENCE: %(prog)s",
                                     add_help=False,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                     )
    # Group 1: General/Default Arguments
    group_default = parser.add_argument_group("General Arguments")
    group_default.add_argument("--gap", type=validate_gap, default=default_gap, metavar="<int>",
                               help="Gap penalty (negative integer)")
    group_default.add_argument("--db", type= check_db_file, default=default_database, metavar="<file>",
                               help="Path to the primary database FASTA file")
    group_default.add_argument("--query", default=default_query, metavar="<file>",
                               help="Path to query FASTA file")
    group_default.add_argument("--method", choices=["blast","sw","stats","tests"], default="blast", metavar="<blast|sw|stats|tests>",
                               help="Choose which search method to run")
    group_default.add_argument("--blosum", choices=[45, 50, 62, 80, 90], type=validate_blosum, default= default_blosum, metavar="<int>",
                               help="Scoring matrix (BLOSUM) to use (45|50|62|80|90)")
    group_default.add_argument("-h", "--help", action="help", help="Show this help message and exit")

    # Group 2: BLAST Arguments
    group_blast = parser.add_argument_group("BLAST-Specific Arguments")
    group_blast.add_argument("--max_alignments", type=validate_max_alignments, default=default_blast_max_alignments, metavar="<int>",
                             help="Maximum number of alignments to output")
    group_blast.add_argument("--max_scores", type=validate_blast_maximum_score, default= default_blast_max_score, metavar="<int>",
                             help= "This determines the maximum number of initial high-scoring hits (HSPs) that the BLAST algorithm internally keeps during the search.")
    group_blast.add_argument("--word_size", type=validate_word_size, default=default_blast_word_size, metavar="<int>",
                             help="Word size for BLAST heuristic")

    # Group 3: Smith-Waterman Arguments
    group_sw = parser.add_argument_group("Smith-Waterman Arguments")
    group_sw.add_argument("--max_sw_scores", type=validate_sw_max_score, default=default_sw_max_score, metavar="<int>",
                          help="This determines the maximum number of initial high-scoring hits (HSPs) that the SW algorithm internally keeps during the search.")

    args = parser.parse_args()

    # Post-parse validations
    if not args.query == default_query:
        args.query = check_query_file(args.query)

    if args.method == "blast":
        if args.max_scores < args.max_alignments:
            print("Error: --max_scores must be >= --max_alignments for BLAST.")
            sys.exit(1)


    return args


def current_arguments():
    args = parse_args()
    print("\nRunning with arguments:")
    for key, value in sorted(vars(args).items()):
        print(f"  {key}: {value}")
    print("\nLOADING...")