# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
#
# - FASTA files to be cleaned
#
# - Metadata file. Mandatory columns: Accession ID, Collection Date, Submission date, Virus name (GISAID)
#                                     Accession (NCBI)
#                                     Headers column (existing or from previous pipeline steps) or list of the columns
#                                     containing the data contained in the sequence header (Personal metadata).
#                                     See "header format".
#
# - header format: If your metadata have not GISAID or Ncbi format, you can use them anyway but you must indicate the
#                  format of the headers of your sequences. If your header has the following format:
#
#                       "Strain/Accession_id/Collection_date"
#
#                  the --header_format string MUST be like this:
#
#                       "\"STRAIN COLUMN NAME\", \"/\" (SEPARATOR), \"ACCESSION ID COLUMN NAME\", \"/\", \"COLLECTION DATE COLUMN NAME\"
#
#                  You have to insert, in order, the name of the columns in the metadata where we can retrieve the
#                  information contained in the header, punctuated by the separators used in the headers.
#                  Another example: if your header has the following format:
#
#                  "Strain|Collection_date"
#
#                  the --header_format string MUST be like this:
#
#                  "\"STRAIN COLUMN NAME\", \"|\" (SEPARATOR), \"COLLECTION DATE COLUMN NAME\""
#
#                  Please pay attention at the white spaces. If your headers contain white spaces, please remember
#                  that all the information contained in the header after the first white space
#                  will be considered as part of the sequence description when parsed. It is preferable to remove or substitute
#                  all teh white spaces.
#
# - Out folder to save the summary file and the diagnostics (.tsv) files (optional).
#   Default out folder: current folder.
#
# - This script allows the user to restart from a checkpoint if the parameter --is_a_checkpoint equals "yes"
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# Sequence filtering ON THE BASIS OF THE FOLLOWING REQUESTS:
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
import polars as pl
import shutil
import time


# FUNCTIONS ############################################################################################################

def sequence_cleaner(camp, usr_meta, ambiguous, complete, tolerance, path_out, index_col_name):

    meta = pl.read_csv(usr_meta, separator="\t", ignore_errors=True).fill_null(value="Undefined")
    print("METADATA OPENED!")

    print(camp)

    # Input columns checking ##########
    meta_columns = meta.columns

    # Public metadata
    if not index_col_name:
        if "Accession ID" not in meta_columns and "Accession" not in meta_columns:
            raise ValueError(
                "Wrong column names or missing columns. Minimum mandatory column for GISAID or NCBI metadata:"
                " see info at --help")

        # GISAID metadata
        elif "Accession ID" in meta_columns:
            type = "GISAID"
            print(type)

        # NCBI metadata
        elif "Accession" in list(meta.columns):
            type = "NCBI"
            print(type)

    # Personal metadata
    else:
        type = "personal"
        print(type)

    # FASTA file type verification and headers making ##########
    # GISAID ---------
    if type == "GISAID":

        with open(camp, "r") as temp:
            for line in temp:
                if "EPI" in line:
                    type = "A"
                    print(type)

                else:
                    type = "B"
                    print(type)

                break

        if type == "B" and "Headers" not in meta.columns:

            meta = meta.with_columns((pl.col("Virus name") + "|"
                                    + pl.col("Collection date")+ "|"
                                    + pl.col("Submission date")).alias('Headers'))

            meta = meta.with_columns(pl.col("Headers").str.replace_all("|Headers", ""))
            meta = meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
            meta = meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
            meta = meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))

        elif type == "A" and "Headers" not in meta.columns:

            meta = meta.with_columns((pl.col("Virus name") + "|"
                                    + pl.col("Accession ID")+ "|"
                                    + pl.col("Collection date")).alias('Headers'))

            meta = meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
            meta = meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
            meta = meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))
            meta = meta.with_columns(pl.col("Headers").str.replace_all("|Headers", ""))

            print("Headers ready")

    # Parameters to count ##########
    tot_written = 0
    file_counter = 0
    length_counter = 0
    iupac_counter = 0
    written_headers = []
    # ##############################

    name = camp.split("/")[-1].replace(".fasta", "")
    file_counter += 1

    with open(f"{path_out}CL_squared_step_02/{name}_cleaned_summary.txt", "w") as doc:
        with open(f"{path_out}CL_squared_step_02/Temp_FASTA/{name}_cleaned_step_02.fasta", "a") as cleaned:

            samples = SeqIO.index(camp, "fasta")

            print("FASTA FILE OPENED!")

            tot_seqs = len(samples)
            print(tot_seqs)

            idx = -1

            # FASTA file type verification and headers making ##########
            # GISAID ---------
            if type == "GISAID":

                with open(camp, "r") as temp:
                    for line in temp:
                        if "EPI" in line:
                            type = "A"
                            print(type)

                        else:
                            type = "B"
                            print(type)

                        break

                if type == "B" and "Headers" not in meta.columns:

                    meta = meta.with_columns((pl.col("Virus name") + "|"
                                            + pl.col("Collection date")+ "|"
                                            + pl.col("Submission date")).alias('Headers'))

                    meta = meta.with_columns(pl.col("Headers").str.replace_all("|Headers", ""))
                    meta = meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
                    meta = meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
                    meta = meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))

                elif type == "A" and "Headers" not in meta.columns:

                    meta = meta.with_columns((pl.col("Virus name") + "|"
                                            + pl.col("Accession ID")+ "|"
                                            + pl.col("Collection date")).alias('Headers'))

                    meta = meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
                    meta = meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
                    meta = meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))
                    meta = meta.with_columns(pl.col("Headers").str.replace_all("|Headers", ""))

                print("Headers ready")

            # NCBI ----------
            elif type == "NCBI":

                meta = meta.with_column(pl.col("Accession").str.replace_all(" ", "_"))
                print("Headers ready")

            # Personal ----------
            elif type == "personal":
                if "Headers" not in meta.columns:

                    if not index_column:
                        raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the index column")

                    if not header_format:
                        raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the header format")

                    header_elements = header_format.split(",")

                    inspected = 0

                    all_idxs = list(cleaned_meta.get_column(index_column))

                    all_seqs = []
                    for idx in all_idxs:
                        inspected += 1
                        header = ""
                        for el in header_elements:
                            if el not in cleaned_meta.columns:
                                header += el
                            else:
                                my_idx = all_idxs.index(idx)
                                temp_col = cleaned_meta.get_column(el)
                                header += f"{temp_col[my_idx]}"

                        all_seqs.append(header)

                    meta = meta.with_columns(pl.Series(name="Headers", values=predictions))


            # SEQUENCES PARSING AND CLEANING ##########
            headers = meta.get_column("Headers")

            for header in headers:
                try:
                    record = samples[header]

                    idx += 1

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
        doc.write(f"Tot files parsed: {file_counter}\n")
        doc.write(f"Tot incomplete sequences: {length_counter}\n")
        doc.write(f"Tot sequences containing non-IUPAC characters: {iupac_counter}\n")
        doc.write(f"Tot written sequences: {tot_written}\n")

        print("FASTA files written")

        # Metadata file ----------

        meta = meta.filter(pl.col("Headers").is_in(written_headers))
        print(f"Metadata length after sequence cleaning: {len(meta)}\n")

        meta.write_csv(f"{path_out}CL_squared_step_02/Temp_META_2/{name}_metadata_cleaned_step_02.tsv", separator='\t')

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

parser = argparse.ArgumentParser(description="Sequence filtering on the basis of the user requests:\n- COMPLETENESS "
                                             "filter: sequences shorter than the required number of bases are discarded"
                                             " (i.e. 29000 bases is the suggested length for SARS-CoV-2 sequences).\n"
                                             "- AMBIGUOUS CHARACTERS filter: sequences with at least one ambiguous "
                                             "character are discarded. The user can set a \"tolerance\" parameter to "
                                             "chose the allowed percentage of ambiguous characters in the sequences.\n",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="User FASTA files path.\nIf you want to inspect multiple FASTA files,"
                                                 " insert the folder that contains all your files.\nATTENTION: IF THE "
                                                 "CHOSEN FOLDER CONTAINS EXTRA FASTA FILES, THEY WILL BE INSPECTED BY "
                                                 "THE SCRIPT")

# User metadata file
parser.add_argument("usr_meta", type=str, help="User metadata file path.\nThe metadata have to be correctly "
                                               "associated with the sequences that you want to inspect.\nIf you have "
                                               "multiple metadata files, concatenate or merge them.\nThe metadata"
                                               " are used in order to quickly access the FASTA file.\n    - If you work"
                                               " with NCBI or GISAID sequences the header of your sequences will be "
                                               "reconstructed starting from the metadata.\n    - If you have personal "
                                               "metadata and sequences format PLEASE ADD A COLUMN CONTAINING THE HEADER"
                                               " OF THE SEQUENCES IN YOUR METADATA FILE AND CALL IT \"Headers\". "
                                               "OTHERWISE, DO PROVIDE THE NAME OF THE COLUMNS CONTAINING THE INFORMATION "
                                               "OF THE HEADER THROUGH THE --header_format PARAMETER.")

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
parser.add_argument("--index_col_name", type=str, help="If your metadata file does not match the GISAID and NCBI "
                                                       "formats and they do not contain the \"Headers\" column, you "
                                                       "must provide the name of the index column and the sequence header"
                                                       "format. See --header_format.")

# Header format for metadata different from GISAID, Ncbi
parser.add_argument("--header_format", type=str, help="If your metadata have not GISAID or Ncbi format and they do not "
                                                      "contain the \"Headers\" column, you can use them anyway but you "
                                                      "must indicate the index column and the format of the headers of "
                                                      "your sequences.\nIf your header has the following format:\n\n\">"
                                                      "Strain/Accession_id/Collection_date\"\n\nthe \"--header_format\""
                                                      " string MUST be like this:\n\n\"\"STRAIN COLUMN NAME\", "
                                                      "\"/\" (SEPARATOR), \"ACCESSION ID COLUMN NAME\", \"/\" "
                                                      "(SEPARATOR), \"COLLECTION DATE COLUMN NAME\"\".\n\nYou basically"
                                                      " have to insert, in order, the names of the columns in the "
                                                      "metadata where we can retrieve the information contained in the "
                                                      "header punctuated by the separators used in the headers.\nAnother"
                                                      " example: if your header has the following format:\n\n \">Strain|"
                                                      "Collection_date\"\n\nthe \"--header_format\" string MUST be "
                                                      "like this:\n\n\"\"STRAIN COLUMN NAME\", \"|\" (SEPARATOR), \""
                                                      "COLLECTION DATE COLUMN NAME\"\n\nPlease pay attention at the white"
                                                      " spaces. If your headers contain white spaces, please remember that"
                                                      " all the information contained in the header after the first white "
                                                      "space will be considered as part of the sequence description when "
                                                      "parsed. It is preferable to remove or substitute all teh white spaces.")

args = parser.parse_args()

#  INPUTS CONTROLLER ----------
if args.usr_fastas.endswith(".fasta"):
    camps = [args.usr_fastas]
elif args.usr_fastas.endswith("/"):
    camps = [i for i in glob.glob(f"{args.usr_fastas}*.fasta")]
else:
    camps = [i for i in glob.glob(f"{args.usr_fastas}/*.fasta")]

print(camps)

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

if not os.path.isdir(f"{path_out}CL_squared_step_02/Temp_META_2"):
    os.makedirs(f"{path_out}CL_squared_step_02/Temp_META_2")

processes = []
for camp in camps:
    p = multiprocessing.Process(target=sequence_cleaner, args=(camp,
                                                                   args.usr_meta,
                                                                   args.ambiguous,
                                                                   args.complete,
                                                                   args.tolerance,
                                                                   path_out,
                                                                   args.index_col_name))
    processes.append(p)
    p.start()

for process in processes:
    process.join()

all_mini_meta = [i for i in glob.glob(f"{path_out}CL_squared_step_02/Temp_META_2/*.tsv")]

new_mm = pd.read_csv(all_mini_meta[0], sep="\t")
print(new_mm)
all_meta_seqs = len(new_mm)
for mm in all_mini_meta[1:]:
    meta_file = pd.read_csv(mm, sep="\t")
    all_meta_seqs += len(meta_file)
    new_mm = pd.concat([new_mm, meta_file], sort=False)
print(len(new_mm))

new_mm.to_csv(f"{path_out}CL_squared_step_02/Cleaned_meta_step_02_2.tsv", sep="\t")

shutil.rmtree(f"{path_out}CL_squared_step_02/Temp_META_2")

