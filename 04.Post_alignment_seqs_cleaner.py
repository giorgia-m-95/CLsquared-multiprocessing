# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - User FASTA files, aligned VS the pathogen reference sequence, to be cleaned
#
# - Metadata file. Mandatory columns: Accession ID, Collection Date, Submission date, Virus name (GISAID)
#                                     Accession (NCBI - please download the accession ID WITH VERSION)
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
# - Admitted percentage of ambiguous characters in the genomic sequences (default 0)
#
# - Number of admitted number of gaps in the head and/or in the tail portion of the sequence (default 90)
#
# - Coding regions file: tsv file containing the name and the genomic locations of the coding regions of the pathogen genome
#			 in a list format ([x, y])
#
# - Out folder for the summary file and the cleaned metadata and FASTA files (optional).
#   Default out folder: current folder.
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# - FASTA aligned sequence filtering.
#   If required by the user, the sequences will be discarded on the basis of the following criteria:
#    - A percentage of ambiguous characters
#    - Head/tail gap trails longer than x bases (default: 0)
#    - Gaps causing frameshifting
# - Metadata cleaning (metadata linked to deleted sequences are removed)
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - Cleaned metadata .tsv file
# - Cleaned FASTA files
# - Summary .txt file
# ----------------------------------------------------------------------------------------------------------------------
#
# ######################################################################################################################

import argparse
from argparse import RawTextHelpFormatter
import ast
from Bio import SeqIO
import glob
import polars as pl
import multiprocessing
from multiprocessing import Process, Lock
import os
import pandas as pd
import re
import shutil


# FUNCTIONS ############################################################################################################
# Pathogen coding regions ----------
def portions_dict(pathogen_portions_file):

    pl.Config.set_tbl_formatting("ASCII_FULL_CONDENSED")

    portions_df = pl.read_csv(pathogen_portions_file, separator='\t', ignore_errors=True).fill_null(value="Undefined")

    portions = {}

    portion_names = list(portions_df.get_column("Name"))

    if len(portion_names) != len(list(set(portion_names))):
        raise ValueError("ATTENTION! Duplicate portion names")
    else:
        for pn in portion_names:
            idx = portion_names.index(pn)
            portions[pn] = portions_df.item(idx, "Positions")

    print(portions)
    return portions


def aligned_seqs_cleaner(camp, meta, ambiguous, tolerance, trails, head, tail, frameshifting, portions, path_out, ref_path):
    print("ENTERED ##########")

    name = camp.split("/")[-1].replace(".fasta", "")

    pl.Config.set_tbl_formatting("ASCII_FULL_CONDENSED")

    with open(f"{path_out}CL_squared_step_04/{name}_sequence_cleaning_summary.txt", "w") as summary:
        # Parameters to count ----------
        disc_ambiguous = 0
        disc_head_tail = 0
        disc_out_frame = 0
        epi_ok = []
        in_seqs = 0
        tot_written_seqs = 0
        # ------------------------------
        # SEQUENCES PARSING AND CLEANING ##########
        headers = meta.get_column("Headers")

        if ref_path:
            ref_seq = SeqIO.read(ref_path, "fasta")
            ref_id = ref_seq.id

        all_msa = [camp]

        for msa in all_msa:
            name = msa.split("/")[-1].replace(".fasta", "")

            algn = SeqIO.index(msa, "fasta")

            summary.write("\n{}: {} starting sequences\n".format(msa.split("/")[-1], len(algn)))

            with open(f"{path_out}CL_squared_step_04/{name}_cleaned_step_04.fasta", "w") as doc:
                try:
                    ref = algn[ref_id]
                    SeqIO.write(ref, doc, "fasta")
                except KeyError:
                    print("REFERENCE IS NOT PRESENT IN THE FASTA FILE")

                for header in headers:
                    try:
                        record = algn[header]
                        sec = record.seq
                        in_seqs += 1

                        # If required, sequences containing ambiguous characters are discarded
                        if ambiguous == "yes":
                            # Maximum percentage of allowed ambiguous characters
                            if tolerance:
                                tolerance_thresh = tolerance
                            else:
                                tolerance_thresh = 0

                            # Percentage of ambiguous characters
                            percentage = ambiguous_characters_checker(record)

                            if percentage <= tolerance_thresh:
                                # If requierd, head/tail gap trails check
                                if trails == "yes":
                                    if head:
                                        start = head
                                    else:
                                        start = 91
                                    if tail:
                                        stop = tail
                                    else:
                                        stop = 91
                                    # Sequences containing gap trails longer than x bases are discarded
                                    if sec[:start] == '-' * start or sec[-stop:] == '-' * stop:
                                        disc_head_tail += 1
                                    else:
                                        # If required, gaps causing frameshift check
                                        if frameshifting == "yes":
                                            out_flag = frameshifting_checker(sec, portions)
                                            if out_flag == 1:
                                                disc_out_frame += 1
                                            elif out_flag == 0:
                                                SeqIO.write(record, doc, "fasta")
                                                tot_written_seqs += 1
                                                epi_ok.append(header)
                                        else:
                                            SeqIO.write(record, doc, "fasta")
                                            tot_written_seqs += 1
                                            epi_ok.append(header)
                                else:
                                    if frameshifting == "yes":
                                        out_flag = frameshifting_checker(sec, portions)
                                        if out_flag == 1:
                                            disc_out_frame += 1
                                        elif out_flag == 0:
                                            SeqIO.write(record, doc, "fasta")
                                            tot_written_seqs += 1
                                            epi_ok.append(header)
                                    else:
                                        SeqIO.write(record, doc, "fasta")
                                        tot_written_seqs += 1
                                        epi_ok.append(header)
                            else:
                                disc_ambiguous += 1
                        else:
                            if trails == "yes":
                                if head:
                                    start = head
                                else:
                                    start = 91
                                if tail:
                                    stop = tail
                                else:
                                    stop = 91
                                # Sequences containing gap trails longer than 90 bases are discarded
                                if sec[:start] == '-' * start or sec[-stop:] == '-' * stop:
                                    disc_head_tail += 1
                                else:
                                    if frameshifting == "yes":
                                        out_flag = frameshifting_checker(sec, portions)
                                        if out_flag == 1:
                                            disc_out_frame += 1
                                        elif out_flag == 0:
                                            SeqIO.write(record, doc, "fasta")
                                            tot_written_seqs += 1
                                            epi_ok.append(header)
                                    else:
                                        SeqIO.write(record, doc, "fasta")
                                        tot_written_seqs += 1
                                        epi_ok.append(header)
                            else:
                                if frameshifting == "yes":
                                    out_flag = frameshifting_checker(sec, portions)
                                    if out_flag == 1:
                                        disc_out_frame += 1
                                    elif out_flag == 0:
                                        SeqIO.write(record, doc, "fasta")
                                        tot_written_seqs += 1
                                        epi_ok.append(header)
                                else:
                                    SeqIO.write(record, doc, "fasta")
                                    tot_written_seqs += 1
                                    epi_ok.append(header)

                    except KeyError:
                        print(f"Header {header} not found in FASTA file")

        # Print reports ----------
        summary.write("\nTOT STARTING SEQUENCES: {}\n".format(in_seqs))
        summary.write("TOT WRITTEN SEQUENCES: {}\n".format(tot_written_seqs))
        if trails == "yes":
            summary.write('SEQS WITH HEAD-TAIL GAPS: {}\n'.format(disc_head_tail))
        if frameshifting == "yes":
            summary.write('FRAMESHIFT SEQS: {}\n'.format(disc_out_frame))

        # Removal of metadata corresponding to discarded seqs ----------
        print("# ###########################################################################################")
        metaf = meta.filter(pl.col("Headers").is_in(epi_ok))
        print("Final metadata length: {}\n".format(len(metaf)))

        if len(metaf) != tot_written_seqs:
            raise ValueError("ERROR! Final metadata do not correspond to written sequences")

        cname = camp.split("/")[-1].replace(".fasta", "")
        metaf.write_csv(f"{path_out}CL_squared_step_04/Temp_META_4/{cname}_step_04.tsv", separator="\t")

    return


# Sequences containing gaps causing frameshift are discarded ----------
def frameshifting_checker(sec, portions):

    flag = 0

    for orf in portions.keys():
        pos_list = ast.literal_eval(portions[orf])
        start = int(pos_list[0]) - 1
        stop = int(pos_list[1])
        temp = sec[start:stop]

        # If the start/stop codon contains at least one gap the sequence is discarded
        if "-" in temp[0:3] or "-" in temp[-3:]:
            flag = 1
            break

        else:
            # Gap groups analysis
            for re.Match in re.finditer(r'[-]+', str(temp)):
                x = list(re.Match.span())
                counter = x[1] - x[0]
                if counter % 3 != 0:
                    flag = 1
                    break
            if flag == 1:
                break
    return flag


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

parser = argparse.ArgumentParser(description="FASTA aligned sequence filtering.\nIf required by the user, the sequences"
                                             " will be discarded on the basis of the following criteria:\n- Presence of"
                                             " ambiguous characters\n- Head/tail gap trails longer than x bases ("
                                             "default: 0)\n- Gaps causing frameshift in the coding regions\nMetadata "
                                             "cleaning: metadata linked to deleted sequences are removed.\nATTENTION: if"
                                             " you want to work in multiprocessing mode but you have a unique FASTA file,"
                                             " please use the 02.1.Headers_to_parallel.py and the 02.1.Parallel_FASTA_maker.sh"
                                             " scripts to create a number of FASTA files corresponding to the number of "
                                             "threads that you have available.",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="User FASTA file's path.\nIf you want to inspect multiple FASTA files,"
                                                 " write the name of the folder that contains all your files.\nEvery "
                                                 "sequence has to be associated to the corresponding metadata "
                                                 "(usr_meta).\nATTENTION: IF THE CHOSEN FOLDER CONTAINS EXTRA FASTA "
                                                 "FILES, THEY WILL BE INSPECTED BY THE SCRIPT.")

# User metadata file
parser.add_argument("usr_meta", type=str, help="User metadata file's path..\nIf you have "
                                               "multiple metadata files, concatenate or merge them.\nIf you work with "
                                               "NCBI meta please download the accession ID WITH VERSION\nIf you work"
                                               " with NCBI or GISAID sequences the header of your sequences will be "
                                               "reconstructed starting from the metadata.\n    - If you have personal "
                                               "metadata and sequences format PLEASE ADD A COLUMN CONTAINING THE HEADER"
                                               " OF THE SEQUENCES IN YOUR METADATA FILE AND CALL IT \"Headers\". "
                                               "OTHERWISE, DO PROVIDE THE NAME OF THE COLUMNS CONTAINING THE INFORMATION "
                                               "OF THE HEADER THROUGH THE --header_format PARAMETER.")

# Ambiguous character check
parser.add_argument("ambiguous", type=str, help="Sequences are inspected for the presence of ambiguous characters.\n"
                                                "If you chose this parameter to be \"yes\" and you do not set the "
                                                "\"tolerance\" parameter, all the sequences with at least one ambiguous"
                                                " character will be discarded.\nIf you want to set a percentage of "
                                                "allowed ambiguos characters set the \"tolerance\" parameter.\nIf you"
                                                " do not want to filter your sequences on the basis of their ambiguous "
                                                "characters, set this parameter as \"no\"")

# Tolerance threshold for ambiguous character checking
parser.add_argument("--tolerance", type=float, help="Percentage of genomic positions allowed to be ambiguous characters"
                                                    " (optional).\nMandatory format: int or float type. The percentage "
                                                    "must be comprised in the range [0, 100]. Default value: 0")

# Head/Tail gap trails check
parser.add_argument("trails", type=str, help="Sequences are inspected for the presence of long head/tail gap trails.\n"
                                             "If you want to set a maximum number of allowed gaps use \"head\" and/or"
                                             " \"tail\" parameters.\n")

# Gaps admitted in the head of the sequence
parser.add_argument("--head", type=int, help="Maximum number of gaps allowed in the head of the sequence (optional).\n"
                                             "Default value: 90.")

# Gaps admitted in the tail of the sequence
parser.add_argument("--tail", type=int, help="Maximum number of gaps allowed in the tail of the sequence (optional).\n"
                                             "Default value: 90.")

# Gaps causing frameshift check
parser.add_argument("frameshifting", type=str, help="Sequences are inspected for the presence of deletions causing "
                                                    "frameshifts in the coding regions.\nAllowed values: yes/no")

# Pathogen coding regions
parser.add_argument("pathogen_portions_file", type=str, help="Path of the .csv/.tsv file containing the name of the "
                                                             "pathogen genome coding regions, each one associated "
                                                             "with the corresponding genomic positions. The file MUST "
                                                             "contain two columns: the first one named \"Name\" "
                                                             "containing the name of each coding region, the second one"
                                                             " named \"Positions\" containing the GENOMIC positions of"
                                                             " each coding region in a list format ([x, y]).")

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

# Reference sequence path
parser.add_argument("--ref_path", type=str, help="Reference sequence path (optional).\nUse this parameter if "
                                                 "your FASTA files contain the reference sequence.")

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

args = parser.parse_args()

if __name__ == "__main__":

    #  INPUTS CONTROLLER ----------
    outfolder = args.outfolder
    if outfolder:
        if outfolder.endswith("/"):
            path_out = outfolder
        else:
            path_out = "{}/".format(outfolder)
    else:
        path_out = ""

    usr_fastas = args.usr_fastas
    if usr_fastas.endswith(".fasta"):
        all_msa = [usr_fastas]
    elif usr_fastas.endswith("/"):
        all_msa = [i for i in glob.glob("{}*.fasta".format(usr_fastas))]
    else:
        all_msa = [i for i in glob.glob("{}/*.fasta".format(usr_fastas))]

    # Metadata controller ----------
    meta = pl.read_csv(args.usr_meta, separator="\t", quote_char=None, ignore_errors=True).fill_null(value="Undefined")

    # Input columns checking ##########
    meta_columns = meta.columns

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

        # NCBI metadata
        elif "Accession" in meta.columns:
            type = "NCBI"
            print(type)

    # Personal metadata
    else:
        type = "personal"
        print(type)

    # Headers making ##########
    index_col_name = args.index_col_name

    # GISAID ---------
    if type == "GISAID":
        camp = all_msa[0]
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
        if "Headers" not in meta.columns:
            meta = meta.with_columns(pl.col("Accession").str.replace_all(" ", "_"))
            print("Headers ready")

    # Personal ----------
    elif type == "personal":
        if "Headers" not in meta.columns:

            if not index_col_name:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the index column")

            if not args.header_format:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the header format")

            header_elements = args.header_format.split(",")

            inspected = 0

            all_idxs = list(meta.get_column(index_col_name))

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

    if not os.path.isdir("{}CL_squared_step_04".format(path_out)):
        os.makedirs("{}CL_squared_step_04".format(path_out))

    if not os.path.isdir(f"{path_out}CL_squared_step_04/Temp_META_4/"):
        os.makedirs(f"{path_out}CL_squared_step_04/Temp_META_4/")

    if args.frameshifting not in ["yes", "no"]:
        raise ValueError("You have entered the wrong input. User-allowed inputs: yes, no")

    if args.ambiguous not in ["yes", "no"]:
        raise ValueError("You have entered the wrong input. User-allowed inputs: yes, no")

    if args.ambiguous == "no" and args.tolerance:
        raise ValueError("You have set \"ambiguous\" parameter as \"no\": set ambiguous parameter as \"yes\" if you "
                         "want to filter your sequences on the basis of the ambiguous characters.")

    if args.trails == "no" and args.head or args.trails == "no" and args.tail:
        raise ValueError("You have set \"trail\" parameter as \"no\": set \"trail\" parameter as \"yes\" if you want to"
                             " filter your sequences on the basis of the head/tail gap trails.")

    portions = portions_dict(args.pathogen_portions_file)

    if os.path.isfile(f"{path_out}CL_squared_step_03/Sequence_cleaning_summary.txt"):
        os.remove(f"{path_out}CL_squared_step_03/Sequence_cleaning_summary.txt")

    # MULTIPROCESSING -------------
    multiprocessing.set_start_method('spawn')

    if len(all_msa) == 1:
        aligned_seqs_cleaner(all_msa[0], meta, args.ambiguous, args.tolerance, args.trails, args.head, args.tail,
                             args.frameshifting, portions, path_out, args.ref_path)

    else:
        lock = Lock()

        processes = []
        for camp in all_msa:

            p = multiprocessing.Process(target=aligned_seqs_cleaner, args=(camp,
                                                                           meta,
                                                                           args.ambiguous,
                                                                           args.tolerance,
                                                                           args.trails,
                                                                           args.head,
                                                                           args.tail,
                                                                           args.frameshifting,
                                                                           portions,
                                                                           path_out,
                                                                           args.ref_path))
            processes.append(p)

        for process in processes:
            process.start()

        for process in processes:
            process.join()

    all_mini_meta = [i for i in glob.glob(f"{path_out}CL_squared_step_04/Temp_META_4/*.tsv")]

    new_mm = pd.read_csv(all_mini_meta[0], sep="\t")
    print(new_mm)
    all_meta_seqs = len(new_mm)
    for mm in all_mini_meta[1:]:
        meta_file = pd.read_csv(mm, sep="\t")
        all_meta_seqs += len(meta_file)
        new_mm = pd.concat([new_mm, meta_file], sort=False)
    print(len(new_mm))

    new_mm.to_csv(f"{path_out}CL_squared_step_04/Cleaned_meta_step_04.tsv", sep="\t")

    shutil.rmtree(f"{path_out}CL_squared_step_04/Temp_META_4")
