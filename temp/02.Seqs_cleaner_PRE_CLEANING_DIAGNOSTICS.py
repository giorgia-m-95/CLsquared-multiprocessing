# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
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
#                  Please pay attention at the white spaces. If your headers contain white spaces, please
#                  remember that all the information contained in the header after the first white space
#                  will be considered as part of the sequence description when parsed. It is preferable to remove or substitute
#                  all teh white spaces.
#
# - Out folder to save the summary file and the diagnostics (.tsv) files (optional).
#   Default out folder: current folder.
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# The script inspects the sequences to filter in order to provide useful information to consciously clean your dataset.
# Specifically, information with respect to the parameters that you will choose in the cleaning step is contained in the
# output (.tsv) file.
# - Sequences length and completeness (sequences characterized by a number of bases >= the required one are considered
# complete)
# - Ambiguous characters: the percentage of ambiguous characters contained in each sequence is provided.
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - Summary (.txt) file
# - Diagnostics (.tsv) files
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

def sequences_inspection(camp, headers, completeness, path_out):

    ffile = camp

    name = ffile.split("/")[-1].replace(".fasta", "")

    with open(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS/Sequences_inspection_diagnostics_summary_{name}.txt", "w") as summary:

        # SEQUENCES PARSING AND CLEANING ##########
        summary.write(f"Working on file {ffile}\n")

        # File indexing
        samples = SeqIO.index(ffile, "fasta")

        tot_seqs = 0
        seqs_dict = {}
        for header in headers:
            # Could not be complete concordance between metadata and FASTA files
            try:
                record = samples[header]

                tot_seqs += 1

                seqs_dict[record.id] = {}

                # Length and completeness check
                seqs_dict[record.id]["Length"] = len(record.seq)

                if len(record.seq) >= completeness:
                    seqs_dict[record.id]["Completeness"] = "yes"
                else:
                    seqs_dict[record.id]["Completeness"] = "no"

                # No IUPAC characters counter
                ambiguous_ch = list(set(list(record.seq)) - set(["A", "T", "C", "G"]))

                all_amb = 0
                for amb in ambiguous_ch:
                    str_amb = record.seq.count(amb)
                    all_amb += str_amb
                percentage = round(all_amb / len(record.seq) * 100, 3)
                seqs_dict[record.id]["Ambiguous characters perc."] = f"{percentage}%"

            except KeyError:
                continue

        summary.write(f"TOT INSPECTED SEQUENCES: {tot_seqs}")

        df = pd.DataFrame.from_dict(seqs_dict).T
        df.to_csv(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS/Temp_META_2/{ffile.split('/')[-1].replace('.fasta', '')}_DIAGNOSTICS.tsv", sep="\t")


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="", formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="Path of the FASTA file that contains the sequences to be "
                                                 "cleaned in CLsquared step 02.\nIf you want to inspect multiple "
                                                 "FASTA files, do insert the folder that contains all your files.\n"
                                                 "ATTENTION: IF THE CHOSEN FOLDER CONTAINS EXTRA FASTA FILES, THEY WILL"
                                                 " BE INSPECTED BY THE SCRIPT.")

# User metadata file
parser.add_argument("usr_meta", type=str, help="User metadata file path.\nThe metadata have to be correctly "
                                               "associated with the sequences that you want to inspect.\nIf you have "
                                               "multiple metadata files, please concatenate or merge them.\nThe metadata"
                                               " are used in order to quickly access the FASTA file.\n    - If you work"
                                               " with NCBI or GISAID sequences the header of your sequences will be "
                                               "reconstructed starting from the metadata.\n    - If you have personal "
                                               "metadata and sequences format PLEASE ADD A COLUMN CONTAINING THE HEADER"
                                               " OF THE SEQUENCES IN YOUR METADATA FILE AND CALL IT \"Headers\". "
                                               "OTHERWISE, DO PROVIDE THE NAME OF THE COLUMNS CONTAINING THE INFORMATION "
                                               "OF THE HEADER THROUGH THE --header_format PARAMETER.")

# Completeness parameter
parser.add_argument("completeness", type=int, help="Number of bases that the genomic sequence is required "
                                                   "to be characterized by in order to be considered complete "
                                                   "(i.e. 29000 bases suggested in the specific case of SARS-CoV-2).")

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
                                                      " spaces. If your headers contain white spaces, please remember that "
                                                      "all the information contained in the header after the first white "
                                                      "space will be considered as part of the sequence description when "
                                                      "parsed. It is preferable to remove or substitute all teh white spaces.")

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: "
                                                  "current folder.")

args = parser.parse_args()


# MAIN #################################################################################################################

#  INPUTS CONTROLLER ----------
# FASTA files
if args.usr_fastas.endswith(".fasta"):
    camps = [args.usr_fastas]
elif args.usr_fastas.endswith("/"):
    camps = [i for i in glob.glob(f"{args.usr_fastas}*.fasta")]
else:
    camps = [i for i in glob.glob(f"{args.usr_fastas}/*.fasta")]

# Outfolder
if args.outfolder:
    outfolder = args.outfolder
    if outfolder.endswith("/"):
        path_out = outfolder
    else:
        path_out = f"{outfolder}/"
else:
    path_out = ""

# Output folder creation ----------
if not os.path.isdir(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS"):
    os.makedirs(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS")

if not os.path.isdir(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS/Temp_META_2"):
    os.makedirs(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS/Temp_META_2")

# Pre-multiprocessing steps ----------
meta = pl.read_csv(args.usr_meta, separator="\t", ignore_errors=True).fill_null(value="Undefined")

meta_columns = meta.columns

if not args.header_format:

    # GISAID metadata
    if "Accession ID" in meta_columns:
        type = "GISAID"
        print(type)

    # NCBI metadata
    elif "Accession" in meta_columns:
        type = "NCBI"
        print(type)

        headers = list(meta.get_column("Accession"))
        print("Headers ready")

    # Personal metadata
    else:

        type = "personal"
        print(type)

        if "Headers" in meta_columns:
            headers = list(meta.get_column("Headers"))
            print("Headers ready")

        else:
            raise ValueError("If your metadata have not GISAID or Ncbi format and they not contain the \"Headers\""
                             "column you must indicate the header format")
else:

    type = "personal"

    if not args.index_col_name:
        raise ValueError("If your metadata have not GISAID or Ncbi format and·they·not·contain·the·\"Headers\""
                         "column you must indicate the index column")
    else:
        header_elements = header_format.split(",")

        inspected = 0

        index_column = args.index_col_name
        all_idxs = list(meta.get_column(index_column))

        headers = []
        for idx in all_idxs:
            inspected += 1
            header = ""
            for el in header_elements:
                if el not in meta.columns:
                    header += el
                else:
                    my_idx = all_idxs.index(idx)
                    temp_col = meta.get_column(el)
                    header += f"{temp_col[my_idx]}"
            headers.append(header)

        print("Headers ready")

# FASTA file type verification and headers making ###################
# GISAID ---------
if type == "GISAID":
    with open(camps[0], "r") as temp:
        for line in temp:
            if "EPI" in line:
                type = "A"
            else:
                type = "B"
            break

    if type == "A" and "Headers" not in meta.columns:

        meta = meta.with_columns((pl.col("Virus name") + "|"
                                + pl.col("Accession ID")+ "|"
                                + pl.col("Collection date")).alias('Headers'))

        meta = meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
        meta = meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
        meta = meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))
        meta = meta.with_columns(pl.col("Headers").str.replace_all("|Headers", ""))

    elif type == "B" and "Headers" not in meta.columns:

        meta = meta.with_columns((pl.col("Virus name") + "|"
                                + pl.col("Collection date")+ "|"
                                + pl.col("Submission date")).alias('Headers'))

        meta = meta.with_columns(pl.col("Headers").str.replace_all("|Headers", ""))
        meta = meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
        meta = meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
        meta = meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))

    headers = meta.get_column("Headers")
    headers = list(headers)
    print("Headers ready")

meta.clear()
print("Metadata CLOSED")

# MULTIPROCESSING ####################
processes = []
for camp in camps:
    p = multiprocessing.Process(target=sequences_inspection, args=(camp,
                                                                   headers,
                                                                   args.completeness,
                                                                   path_out))
    processes.append(p)
    p.start()

for process in processes:
    process.join()

all_mini_meta = [i for i in glob.glob(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS/Temp_META_2/*.tsv")]

new_mm = pl.read_csv(all_mini_meta[0], separator="\t", ignore_errors=True)
print(new_mm)
all_meta_seqs = len(new_mm)
for mm in all_mini_meta[1:]:
    meta_file = pl.read_csv(mm, separator="\t", ignore_errors=True)
    all_meta_seqs += len(meta_file)
    new_mm = pl.concat([new_mm, meta_file])
print(len(new_mm))

new_mm.write_csv(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS/Cleaned_meta_step_02_2.tsv", separator="\t")

shutil.rmtree(f"{path_out}CL_squared_PRE_step_02-DIAGNOSTICS/Temp_META_2")
