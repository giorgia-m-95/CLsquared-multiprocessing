# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - FASTA files to be grouped
#
# - Metadata file. Mandatory columns: Accession ID, Collection Date, Submission date, Virus name (GISAID)
#                                     Accession (NCBI - please download the accession ID WITH VERSION)
#                                     Headers column (existing or from previous pipeline steps) or list of the columns
#                                     containing the data contained in the sequence header (Personal metadata).
#                                     See "header format".
#
# - Name of the metadata column that contains the information necessary to create the FASTA files containing the grouped
#   sequences (e.g viral variant column name).
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
# - Out folder to save the summary file and the cleaned FASTA files (optional). Default out folder: current folder.
#
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# Creation of FASTA files containing the starting sequences divided in groups on the basis of the viral variant or
# strain. The user indicates the metadata column containing the information for the division criteria. This script is
# recommended if the user wants to align specific sequences of the starting dataset versus a chosen reference sequence
# (e.g. starting dataset containing multiple viral strains associated with distinct reference sequences) or if the user
# wants to perform the clustering step 05 as a cleaning step.
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - FASTA files containing the grouped sequences
# - Summary (.txt) file
# ----------------------------------------------------------------------------------------------------------------------
#
# ######################################################################################################################

import argparse
import time
from argparse import RawTextHelpFormatter
from Bio import AlignIO, SeqIO
import glob
import polars as pl
import multiprocessing
from multiprocessing import Lock
import os
import pandas as pd
import shutil

# FUNCTIONS ############################################################################################################


def seqs_grouper(camp, meta_path, grouping_column, path_out, unique, ref_path):
    print("ENTERED! ###########################################################")

    name = camp.split("/")[-1].replace(".fasta", "")
    print(name)

    meta = pl.read_csv(meta_path, separator="\t", ignore_errors=True).fill_null(value="Undefined")

    pl.Config.set_tbl_formatting("ASCII_FULL_CONDENSED")
    print(meta)

    cat_col = list(meta.get_column(grouping_column))
    cathegories = list(set(cat_col))
    print(cathegories)

    written_seqs = 0
    for cat in cathegories:
        with open(f"{path_out}CL_squared_step_C1/{name}_{cat}_cleaned_summary.txt", "w") as doc:
            doc.write(f"\n{cat} --------------------\n")

            cat_meta = meta.filter(pl.col(grouping_column) == cat)
            doc.write(f"TOT {cat} SEQUENCES FROM META: {len(cat_meta)}\n")

            written_cat = 0

            if unique == "yes":
                if not os.path.isdir(f"{path_out}CL_squared_step_C1/Temp_{cat}/"):
                    os.makedirs(f"{path_out}CL_squared_step_C1/Temp_{cat}/")
                path = f"{path_out}CL_squared_step_C1/Temp_{cat}/{name}_cleaned_step_C1_{cat}.fasta"
            else:
                path = f"{path_out}CL_squared_step_C1/{name}_cleaned_step_C1_{cat}.fasta"

            if ref_path and unique == "no":
                ref_seq = SeqIO.read(ref_path, "fasta")
                ref_id = ref_seq.id

            with open(path, "w") as cleaned:

                samples = SeqIO.index(camp, "fasta")

                print("FASTA FILE OPENED!")

                tot_seqs = len(samples)
                print(tot_seqs)

                # SEQUENCES PARSING AND CLEANING ##########
                headers = cat_meta.get_column("Headers")

                if ref_path and unique == "no":
                    try:
                        ref = samples[ref_id]
                        SeqIO.write(ref, cleaned, "fasta")
                    except KeyError:
                        print("REFERENCE IS NOT PRESENT IN THE FASTA FILE")
                        raise ValueError("Reference path inserted. Reference sequence not present in the FASTA file.")

                for header in headers:
                    try:
                        record = samples[header]
                        SeqIO.write(record, cleaned, "fasta")
                        written_seqs += 1
                        written_cat += 1
                    except KeyError:
                        print(f"Header {header} not found in FASTA file")

                doc.write(f"WRITTEN SEQS: {written_cat}\n")

    return


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="Creation of FASTA files containing the starting sequences divided in "
                                             "groups on the basis of the viral variant or strain.\nThe user must indicate"
                                             " the metadata column containing the information for the division criteria."
                                             "\nThis script is recommended if the user wants to align specific sequences"
                                             " of the starting dataset versus a chosen reference sequence (e.g. starting"
                                             " dataset containing multiple viral strains associated with distinct reference"
                                             " sequences) or if the user wants to perform the clustering step 05 as a "
                                             "cleaning step.",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="User FASTA files path.\nIf you want to inspect multiple FASTA files,"
                                                 " insert the folder that contains all your files.\nATTENTION: IF THE "
                                                 "CHOSEN FOLDER CONTAINS EXTRA FASTA FILES, THEY WILL BE INSPECTED BY "
                                                 "THE SCRIPT")

# User metadata file
parser.add_argument("usr_metadata", type=str, help="Metadata file path.\nSequence ID and Variant columns are mandatory.\n"
                                                   "If you work with NCBI metadata please download the accession ID WITH"
                                                   " VERSION.\nIf you have multiple metadata files, concatenate or merge"
                                                   " them. MANDATORY FILE FORMAT: .tsv")

# Grouping infomation column
parser.add_argument("grouping_column", type=str, help="Metadata column on the basis of which you want to group your "
                                                      "sequences in distinct FASTA files.")

parser.add_argument("unique", type=str, help="If you want to create a unique FASTA file containing all the sequences "
                                            "belonging to a specific category coming from all your samples, set this"
                                            " parameter as \"yes\", otherwise, set it as \"no\".")

# Output folder path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

# Index column name (sequence id)
parser.add_argument("--index_col_name", type=str, help="If your metadata file does not match the GISAID and NCBI "
                                                       "formats and they do not contain the \"Headers\" column, you "
                                                       "must provide the name of the index column and the sequence header"
                                                       "format. See --header_format.")

# Metadata columns
parser.add_argument("--meta_columns_list", type=str, help="If your metadata file doesn't match the GISAID and NCBI "
                                                          "formats, you MUST provide a list of columns as an extra input"
                                                          ". The list should follow this pattern:\n[\"sequence id "
                                                          "column name\", \"host column name\", \"place column name\","
                                                          " \"collection date column name\"].\n ATTENTION: if you have "
                                                          "a missing column, maintain the order and use \"\" in its "
                                                          "place. For instance, if you lack the \"host\" column, your "
                                                          "list should look like this: [\"sequence id column name\", "
                                                          "\"\", \"place column name\", \"collection date column name"
                                                          "\"]. It is MANDATORY TO HAVE A SEQUENCE ID COLUMN.")

# Reference sequence path
parser.add_argument("--ref_path", type=str, help="Reference sequence path (optional).\nUse this parameter if "
                                                 "your FASTA files contain the reference sequence.")

args = parser.parse_args()


# INPUTS CONTROLLER ----------

if __name__ == "__main__":

    # FASTA FILES ----------
    if args.usr_fastas.endswith(".fasta"):
        camps = [args.usr_fastas]
    elif args.usr_fastas.endswith("/"):
        camps = [i for i in glob.glob(f"{args.usr_fastas}*.fasta")]
    else:
        camps = [i for i in glob.glob(f"{args.usr_fastas}/*.fasta")]
    print(camps)

    # METADATA ----------

    meta = pl.read_csv(args.usr_metadata, separator="\t", ignore_errors=True).fill_null(value="Undefined")
    print("METADATA OPENED!")

    # Input columns checking ##########
    meta_columns = meta.columns

    # Public metadata
    if not args.index_col_name:
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

    # GISAID ---------
    if type == "GISAID":

        camp = camps[0]

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

        meta = meta.with_columns(pl.col("Accession").str.replace_all(" ", "_"))

        if "Headers" not in meta.columns:
            all_seqs = list(meta.get_column("Accession"))
            meta = meta.with_columns(pl.Series(name="Headers", values=all_seqs))

            print("Headers ready")

    # Personal ----------
    elif type == "personal":
        if "Headers" not in meta.columns:

            index_column = args.index_col_name

            if not index_column:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the index column")

            if not args.header_format:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the header format")

            header_elements = args.header_format.split(",")

            inspected = 0

            all_idxs = list(meta.get_column(index_column))

            all_seqs = []
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
                all_seqs.append(header)

            meta = meta.with_columns(pl.Series(name="Headers", values=all_seqs))

    # OUTFOLDER ----------
    if args.outfolder:
        outfolder = args.outfolder
        if outfolder.endswith("/"):
            path_out = outfolder
        else:
            path_out = f"{outfolder}/"
    else:
        path_out = ""

    if not os.path.isdir(f"{path_out}CL_squared_step_C1"):
        os.makedirs(f"{path_out}CL_squared_step_C1")

    meta.write_csv(f"{path_out}CL_squared_step_C1/Cleaned_meta_step_C.tsv", separator="\t")
    meta_path = f"{path_out}CL_squared_step_C1/Cleaned_meta_step_C.tsv"

    # UNIQUE FASTA ----------
    if args.unique not in ["yes", "no"]:
        raise ValueError("Wrong value inserted. Only \"yes\" and \"no\" values are permitted")
    else:
        unique = args.unique
        print(unique)

    # Multiprocessing ----------
    grouping_column = args.grouping_column
    print(grouping_column)

    cat_col = list(meta.get_column(args.grouping_column))
    cathegories = list(set(cat_col))

    multiprocessing.set_start_method('spawn', force=True)

    if len(camps) == 1:
        seqs_grouper(camps[0], meta_path, grouping_column, path_out, unique, args.ref_path)
    else:
        processes = []
        for camp in camps:
            p = multiprocessing.Process(target=seqs_grouper, args=(camp,
                                                                   meta_path,
                                                                   grouping_column,
                                                                   path_out,
                                                                   unique,
                                                                   args.ref_path))
            processes.append(p)

        for process in processes:
            process.start()

        for process in processes:
            process.join()

    if unique == "yes":

        print("WRITING UNIQUE FASTA FILES")

        for cat in cathegories:
            print(f"Working on {cat} ----------")

            all_mini_fasta = [i for i in glob.glob(f"{path_out}CL_squared_step_C1/Temp_{cat}/*.fasta")]
            print(all_mini_fasta)

            with open(f"{path_out}CL_squared_step_C1/{cat}_sequences.fasta", "w") as ffile:

                if args.ref_path:
                    ref_seq = SeqIO.read(args.ref_path, "fasta")
                    SeqIO.write(ref_seq, ffile, "fasta")
                    # ffile.write(f">{ref_seq.id}\n")
                    # ffile.write(f"{ref_seq.seq}\n")

                for amf in all_mini_fasta:
                    seqs = list(SeqIO.parse(amf, "fasta"))
                    SeqIO.write(seqs, ffile, "fasta")
                    # with open(amf, "r") as doc:
                    #     for line in doc:
                    #         ffile.write(line)

            shutil.rmtree(f"{path_out}CL_squared_step_C1/Temp_{cat}")
