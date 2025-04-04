# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - FASTA files to be cleaned.
#
# - Metadata.
#   MANDATORY MINIMUM FORMAT:
#    - Columns: Accession ID (sequence ID) (GISAID metadata format)
#               Accession (sequence ID) (NCBI metadata format)
#               If your metadata file doesn't match the GISAID and NCBI formats you can provide the name of the column
#               containing the SEQUENCE ID THAT IS THE ONLY MANDATORY ONE.
#
# - Out folder to save summary files, cleaned metadata and FASTA files (optional).
#   Default out folder: current folder.
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# Sequence filtering ON THE BASIS OF THE USER REQUESTS:
# - COMPLETENESS filter: sequences shorter than the required number of bases can be discarded (i.e. 29000 bases is
#   the suggested length for SARS-CoV-2 sequences).
# - AMBIGUOUS CHARACTERS filter: sequences with at least one ambiguous character can be discarded. The user can
#   OPTIONALLY set a tolerance parameter to choose the allowed percentage of ambiguous characters in the sequences.
#
# Metadata cleaning: the metadata linked to the deleted sequences are removed.
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - Cleaned metadata (.tsv) file
# - Cleaned FASTA files
# - Summary (.txt) file
# ----------------------------------------------------------------------------------------------------------------------
#
# ######################################################################################################################


import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import glob
import multiprocessing
import os
import pandas as pd
import polars as pl
import shutil


# FUNCTIONS ############################################################################################################

# TODO: ADD THE POSSIBILITY TO FURNISH THE COLUMNS LIST AND THE HEADER FORMAT FOR PERSONAL METADATA --> IMPORT A
#  FUNCTION (?)

def sequence_cleaner(camp, headers, ambiguous, complete, tolerance, path_out):

    # Parameters to count ##########
    tot_written = 0
    length_counter = 0
    iupac_counter = 0
    written_headers = []
    # ##############################

    file = camp

    name = file.split("/")[-1].replace(".fasta", "")

    with open(f"{path_out}CL_squared_step_02/{name}_cleaned_summary.txt", "w") as doc:
        with open(f"{path_out}CL_squared_step_02/Temp_FASTA/{name}_cleaned_step_02.fasta", "a") as cleaned:

            samples = SeqIO.index(file, "fasta")

            print("FASTA FILE OPENED!")

            tot_seqs = len(samples)
            print(tot_seqs)

            # SEQUENCES PARSING AND CLEANING ##########
            for header in headers:

                try:
                    record = samples[header]

                    # If required, completeness check:
                    if complete:
                        if len(record.seq) >= complete:
                            # If required, IUPAC characters check: removal of sequences with characters
                            # different from A,T,C and G
                            if ambiguous == "yes":
                                # Maximum percentage of allowed ambiguous characters
                                if tolerance:
                                    tolerance_thresh = tolerance
                                else:
                                    tolerance_thresh = 0

                                # Percentage of ambiguous characters
                                percentage = ambiguous_characters_checker(record)

                                if percentage <= tolerance_thresh:
                                    # Control OK
                                    SeqIO.write(record, cleaned, "fasta")
                                    written_headers.append(header)
                                    tot_written += 1
                                else:
                                    iupac_counter += 1

                            # In this case sequences will not be checked for ambiguous characters
                            elif ambiguous == "no":
                                SeqIO.write(record, cleaned, "fasta")
                                written_headers.append(header)
                                tot_written += 1
                        else:
                            length_counter += 1

                    else:
                        if ambiguous == "yes":
                            # Maximum percentage of allowed ambiguous characters
                            if tolerance:
                                tolerance_thresh = tolerance
                            else:
                                tolerance_thresh = 0

                            # Percentage of ambiguous characters
                            percentage = ambiguous_characters_checker(record)

                            if percentage <= tolerance_thresh:
                                # Control OK
                                SeqIO.write(record, cleaned, "fasta")
                                written_headers.append(header)
                                tot_written += 1
                            else:
                                iupac_counter += 1

                        # In this case sequences will not be checked for ambiguous characters
                        elif ambiguous == "no":
                            SeqIO.write(record, cleaned, "fasta")
                            written_headers.append(header)
                            tot_written += 1
                except KeyError:
                    print(f"Header {header} not found in FASTA file")

        doc.write(f"Tot sequences: {tot_seqs}\n")
        doc.write(f"Tot incomplete sequences: {length_counter}\n")
        doc.write(f"Tot sequences containing non-IUPAC characters: {iupac_counter}\n")
        doc.write(f"Tot written sequences: {tot_written}\n")

        print("FASTA files written")

        # Metadata file ----------
        with open(f"{path_out}CL_squared_step_02/Logs/{name}_cleaned_summary.log", "w") as log:
            for head in written_headers:
                log.write(f"{head}\n")


# Sequence check for characters ----------
def ambiguous_characters_checker(record):
    ambiguous_ch = list(set(list(record.seq)) - set(["A", "T", "C", "G", "-"]))
    all_amb = 0
    for amb in ambiguous_ch:
        str_amb = record.seq.count(amb)
        all_amb += str_amb
    percentage = all_amb / len(record.seq) * 100

    return percentage


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="Sequence filtering on the basis of the user requests:\n- COMPLETENESS "
                                             "filter: sequences shorter than the required number of bases are discarded"
                                             " (i.e. 29000 bases is the suggested length for SARS-CoV-2 sequences).\n"
                                             "- AMBIGUOUS CHARACTERS filter: sequences with at least one ambiguous "
                                             "character are discarded. The user can set a \"tolerance\" parameter to "
                                             "chose the allowed percentage of ambiguous characters in the sequences.\n"
                                             "Metadata must match the inspected sequences.",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="User FASTA files path.\nIf you want to inspect multiple "
                                                 "FASTA files, insert the folder that contains all your files.\nEvery "
                                                 "sequence has to be associated with the corresponding metadata "
                                                 "(usr_meta).\nATTENTION: IF THE CHOSEN FOLDER CONTAINS EXTRA FASTA "
                                                 "FILES, THEY WILL BE INSPECTED BY THE SCRIPT.")

# User metadata file
parser.add_argument("usr_meta", type=str, help="User metadata file path.\nThe metadata have to be correctly "
                                               "associated with the sequences that you want to inspect.\nIf you have "
                                               "multiple metadata files, concatenate or merge them.")

# Ambiguous characters check
parser.add_argument("ambiguous", type=str, help="Sequences are inspected for the presence of ambiguous characters.\n"
                                                "If you chose this parameter to be \"yes\" and you do not set the "
                                                "\"tolerance\" parameter, all the sequences with at least one ambiguous"
                                                " character will be discarded.\nIf you want to set the percentage of "
                                                "allowed ambiguos characters set the \"tolerance\" parameter.\nIf you "
                                                "do not want to filter your sequences on the basis of their ambiguous "
                                                "characters, set this parameter as \"no\"")

# Completeness check
parser.add_argument("--complete", type=int, help="Sequences completeness check (optional).\nIf you set this parameter, "
                                                 "all the sequences shorter than the number of set bases will be "
                                                 "discarded.")

# Percentage of admitted ambiguous characters
parser.add_argument("--tolerance", type=float, help="Percentage of genomic positions allowed to be ambiguous characters"
                                                    " (optional).\nMandatory format: int or float type. The percentage "
                                                    "must be comprised in the range [0, 100]. Default value: 0")

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")


# Index column name (sequence id)
parser.add_argument("--index_col_name", type=str, help="If your metadata file doesn't match the GISAID and NCBI "
                                                       "formats you can provide the name of the column containing the "
                                                       "SEQUENCE ID THAT IS THE ONLY MANDATORY ONE. If it is unnamed, "
                                                       "use its index number (e.g., 0, 1, ...).")

args = parser.parse_args()

#  INPUTS CONTROLLER ----------
if args.usr_fastas.endswith(".fasta"):
    camps = [args.usr_fastas]
elif args.usr_fastas.endswith("/"):
    camps = [i for i in glob.glob(f"{args.usr_fastas}*.fasta")]
else:
    camps = [i for i in glob.glob(f"{args.usr_fastas}/*.fasta")]

if args.outfolder:
    outfolder = args.outfolder
    if outfolder.endswith("/"):
        path_out = outfolder
    else:
        path_out = f"{outfolder}/"
else:
    path_out = ""

if args.ambiguous not in ["yes", "no"]:
    raise ValueError("Wrong input. User-allowed inputs: yes, no")

if args.ambiguous == "no" and args.tolerance:
    raise ValueError("You set a tolerance parameter: set ambiguous parameter as \"yes\" if you want to filter"
                     " your sequences on the basis of the ambiguous characters.")

# Output folder creation ----------
if not os.path.isdir(f"{path_out}CL_squared_step_02"):
    os.makedirs(f"{path_out}CL_squared_step_02")

if not os.path.isdir(f"{path_out}CL_squared_step_02/Temp_FASTA"):
    os.makedirs(f"{path_out}CL_squared_step_02/Temp_FASTA")

if not os.path.isdir(f"{path_out}CL_squared_step_02/Logs"):
    os.makedirs(f"{path_out}CL_squared_step_02/Logs")

meta = pl.read_csv(args.usr_meta, separator="\t", ignore_errors=True).fill_null(value="Undefined")
print("METADATA OPENED!")

# Input columns checking ##########
meta_columns = meta.columns
print(meta_columns)

# Public metadata
if not args.index_col_name:
    if "Accession ID" not in meta_columns and "Accession" not in meta_columns:
        raise ValueError(
            "Wrong column names or missing columns. Minimum mandatory column for GISAID or NCBI metadata:"
            " Sequence accession ID")

    # GISAID metadata
    elif "Accession ID" in meta_columns:
        type = "GISAID"
        print(type)

    elif "Accession" in list(meta.columns):
        type = "NCBI"
        print(type)
        headers = list(meta.select(pl.col(["Accession"])))
        print("Headers ready")

# Personal metadata
else:
    type = "personal"
    print(type)
    if "Headers" not in meta_columns:
        headers = list(meta.select(pl.col([args.index_col_name])))
        print("Headers ready")
        print(headers)
    else:
        headers = list(meta.select(pl.col(["Headers"])))
        print("Headers ready")

# FASTA file type verification and headers making ##########
# GISAID ---------
if type == "GISAID":
    with open(camps[0], "r") as temp:
        for line in temp:
            if "EPI" in line:
                type = "A"
            else:
                type = "B"
            break

    if type == "A":
        meta = meta.with_columns(pl.concat_str([pl.col("Virus name"), pl.col("Accession ID"),
                                                pl.col("Collection date")], separator="|").alias("Headers"))
        meta.with_columns(pl.col("Headers").str.replace("|Headers", ""))

    elif type == "B":
        meta = meta.with_columns(pl.concat_str([pl.col("Virus name"), pl.col("Collection date"),
                                                pl.col("Submission date")], separator="|").alias("Headers"))
        meta.with_columns(pl.col("Headers").str.replace("|Headers", ""))

    headers = meta.get_column("Headers")
    headers = list(headers)
    print("Headers ready")

# meta.clear()
# print("Metadata CLOSED")

processes = []
for camp in camps:
    p = multiprocessing.Process(target=sequence_cleaner, args=(camp,
                                                               headers,
                                                               args.ambiguous,
                                                               args.complete,
                                                               args.tolerance,
                                                               path_out))
    processes.append(p)
    p.start()

for process in processes:
    process.join()

all_mini_heads = [i for i in glob.glob(f"{path_out}CL_squared_step_02/Logs/*.log")]

all_wr_headers = []
for mm in all_mini_heads:
    with open(mm, "r") as doc:
        for line in doc:
            all_wr_headers.append(line.strip())

meta.row(by_predicate=(pl.col("Headers").isin(all_wr_headers)))

meta.write_csv("{path_out}CL_squared_step_02/Cleaned_metadata.tsv", separator="\t")

shutil.rmtree(f"{path_out}CL_squared_step_02/Logs")
