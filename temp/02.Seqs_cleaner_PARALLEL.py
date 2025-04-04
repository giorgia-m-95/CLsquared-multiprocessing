# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - User's FASTA files to be cleaned:
#   MANDATORY HEADERS FORMAT: >Place/Accession ID (sequence ID)/Collection date
#
#   If your headers are not formatted like this you can input a file with the formatted headers that can be used.
#   The file MUST BE FORMATTED AS FOLLOWS:
#   current_header\tcorresponding_formatted_header
#
#   A script to create this file starting from your metadata and current headers is furnished: "Headers_formatter.py"
#
#   You can alternatively edit all the headers of the FASTA file through the "Headers_creator.py" script.
#   This option is not recommended for large FASTA files.
#
# - User's metadata.
#   MANDATORY MINIMUM FORMAT:
#    - Columns: Accession ID (sequence ID) (GISAID metadata format)
#               Accession (sequence ID) (NCBI metadata format)
#               If your metadata file doesn't match the GISAID and NCBI formats you can provide the name of the column
#               containing the SEQUENCE ID THAT IS THE ONLY MANDATORY ONE. If it is unnamed, use its index number
#               (e.g., 0, 1, ...).
#
# - Out folder to save the summary file and the cleaned metadata and FASTA files (optional).
#   Default out folder: current folder.
#
# - This script allows the user to restart from a checkpoint if the parameter --is_a_checkpoint equals "yes"
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# Sequence filtering ON THE BASIS OF THE USER'S REQUESTS:
# - COMPLETENESS filter: sequences shorter than the required number of bases can be discarded (i.e. 29000 bases is
#   the suggested length for SARS-CoV-2 sequences).
# - AMBIGUOUS CHARACTERS filter: sequences with at least one ambiguous character can be discarded. The user can
#   OPTIONALLY set a tolerance parameter to chose the allowed percentage of ambiguous characters in the sequences.
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
import shutil
import time


# FUNCTIONS ############################################################################################################

def sequence_cleaner(camp, usr_meta, ambiguous, complete, tolerance, path_out, formatted_headers, index_col_name):

    # Sequence check ----------
    meta = pd.read_csv(usr_meta, sep=None, engine="python")
    print("METADATA OPENED!")

    # Input columns checking ----------
    meta_columns = meta.columns
    if not index_col_name:
        if "Accession ID" not in meta_columns and "Accession" not in meta_columns:
            raise ValueError(
                "Wrong column names or missing columns. Minimum mandatory column for GISAID or NCBI metadata:"
                " Sequence accession ID")
        elif "Accession ID" in meta_columns:
            meta.set_index("Accession ID", inplace=True)
        elif "Accession" in list(meta.columns):
            meta.set_index("Accession", inplace=True)
    else:
        if index_col_name.isdigit():
            meta.set_index(int(index_col_name), inplace=True)
        else:
            meta.set_index(index_col_name, inplace=True)

    all_ids = list(meta.index)

    # Parameters to count
    tot_written = 0
    written_epis = []
    file_counter = 0
    length_counter = 0
    iupac_counter = 0
    tot_seqs = 0

    file = camp

    name = file.split("/")[-1].replace(".fasta", "")
    file_counter += 1

    with open(f"{path_out}CL_squared_step_02/{name}_cleaned_summary.txt", "w") as doc:
        with open(f"{path_out}CL_squared_step_02/Temp_FASTA/{name}_cleaned_step_02.fasta", "a") as cleaned:
            samples = list(SeqIO.parse(file, "fasta"))
            print("FASTA FILE OPENED!")

            tot_seqs += len(samples)

            idx = -1

            for record in samples:

                idx += 1

                if not formatted_headers:
                    epi = record.id.split('/')[1]
                else:
                    with open(formatted_headers, "r") as fh:
                        flag = 0
                        for line in fh:
                            if line.split("\t")[0] == record.id:
                                current_header = line.split("\t")[1]
                                record.id = current_header
                                record.description = ""
                                epi = current_header.split('/')[1]
                                flag = 1
                        if flag == 0:
                            print(f"Your current header is not contained in the headers file: {record.id}")
                record.description = ""

                if epi in all_ids:

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
                                    written_epis.append(epi)
                                    tot_written += 1
                                else:
                                    iupac_counter += 1

                            # In this case sequences will not be checked for ambiguous characters
                            elif ambiguous == "no":
                                SeqIO.write(record, cleaned, "fasta")
                                written_epis.append(epi)
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
                                written_epis.append(epi)
                                tot_written += 1
                            else:
                                iupac_counter += 1

                        # In this case sequences will not be checked for ambiguous characters
                        elif ambiguous == "no":
                            SeqIO.write(record, cleaned, "fasta")
                            written_epis.append(epi)
                            tot_written += 1

        doc.write(f"Tot sequences: {tot_seqs}")
        doc.write(f"Tot files parsed: {file_counter}")
        doc.write(f"Tot incomplete sequences: {length_counter}")
        doc.write(f"Tot sequences containing non-IUPAC characters: {iupac_counter}")
        doc.write(f"Tot written sequences: {tot_written}")

        print("FASTA files written")

        # Metadata file ----------
        meta = meta[meta.index.isin(written_epis)]
        print(f"Metadata length after sequence cleaning: {len(meta)}\n")

        meta.to_csv(f"{path_out}CL_squared_step_02/Temp_META/{name}_metadata_cleaned_step_02.tsv", sep='\t')

        # Matching control: all seqs have corresponding metadata and viceversa ----------
        if len(meta) != tot_written:
            raise ValueError("ERROR: metadata length != total number of written sequences")

        time.sleep(2)


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

parser = argparse.ArgumentParser(description="Sequence filtering on the basis of the user's requests:\n- COMPLETENESS "
                                             "filter: sequences shorter than the required number of bases are discarded"
                                             " (i.e. 29000 bases is the suggested length for SARS-CoV-2 sequences).\n"
                                             "- AMBIGUOUS CHARACTES filter: sequences with at least one ambiguous "
                                             "character are discarded. The user can set a \"tolerance\" parameter to "
                                             "chose the allowed percentage of ambiguous characters in the sequences.\n"
                                             "Metadata must match the inspected sequences.",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="User FASTA files path.\nIf you want to inspect multiple FASTA files,"
                                                 " insert the folder that contains all your files.\nEvery "
                                                 "sequence has to be associated with the corresponding metadata "
                                                 "(usr_meta).\nATTENTION: IF THE CHOSEN FOLDER CONTAINS EXTRA FASTA "
                                                 "FILES, THEY WILL BE INSPECTED BY THE SCRIPT.\nMandatory sequence's "
                                                 "header format: >Place/Accession ID (sequence ID)/Collection date.\nIf"
                                                 " the format of your headers does not match this one, a formatting "
                                                 "script is provided:\nHeader_formatter.py")

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

# Formatted headers file
parser.add_argument("--formatted_headers", type=str, help="File containing the formatted headers path.\n"
                                                          "The format of the rows of the file must be as follows:"
                                                          "CURRENT_SEQUENCE_HEADER\tNEW_HEADER.\n"
                                                          "A script to format your headers like this is furnished and "
                                                          "is named 'Headers_formatter.py'")

# Checkpoint
parser.add_argument("--is_a_checkpoint", type=str, help="If the script is restarting from a checkpoint set this \n"
                                                        "parameter equal to 'yes' (optional).")

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

if not os.path.isdir(f"{path_out}CL_squared_step_02/Temp_META"):
    os.makedirs(f"{path_out}CL_squared_step_02/Temp_META")

processes = []
for camp in camps:
    p = multiprocessing.Process(target=sequence_cleaner, args=(camp,
                                                                   args.usr_meta,
                                                                   args.ambiguous,
                                                                   args.complete,
                                                                   args.tolerance,
                                                                   path_out,
                                                                   args.formatted_headers,
                                                                   args.index_col_name))
    processes.append(p)
    p.start()

for process in processes:
    process.join()

all_mini_meta = [i for i in glob.glob(f"{path_out}CL_squared_step_02/Temp_META/*.tsv")]

new_mm = all_mini_meta[0]
for mm in all_mini_meta[1:]:
    new_mm = new_mm.append(mm, ignore_index=True)

new_mm.to_csv(f"{path_out}CL_squared_step_02/Cleaned_meta_step_02.tsv", sep="\t")

shutil.rmtree(f"{path_out}CL_squared_step_02/Temp_META")