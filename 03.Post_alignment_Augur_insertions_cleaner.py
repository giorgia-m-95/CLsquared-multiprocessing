# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - .csv files where insertions have been annotated --> augur align pipeline out files
#
# - User FASTA files, aligned VS the pathogen reference sequence or original FASTA files before alignment VS the pathogen
#   reference sequence, to be cleaned
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
# - Out folder for the summary file and the cleaned metadata (optional).
#   Default out folder: current folder.
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# The sequences are cleaned on the basis of the insertions that were trimmed by the augur align pipeline.
# Sequences with out of frame insertions in coding regions are deleted.
# The user can decide to only filter his sequences or, additionally, to also trim their extremities.
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - FASTA files with the cleaned sequences.
# - Cleaned metadata .tsv file.
# - Summary .txt file.
# ----------------------------------------------------------------------------------------------------------------------
#
# ######################################################################################################################

import argparse
from argparse import RawTextHelpFormatter
import ast
from Bio import SeqIO
import multiprocessing
import glob
import os
import pandas as pd
import polars as pl
import shutil


# FUNCTIONS ############################################################################################################
# Pathogen coding regions ----------
def portions_dict(pathogen_portions_file):
    portions_df = pl.read_csv(pathogen_portions_file, separator="\t", ignore_errors=True).fill_null(value="Undefined")

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


# Insertions inspection ----------
def insertions_inspector(camp, portions, path_out):

    insertion_files_path = camp

    if insertion_files_path.endswith(".fasta"):
        ins_files = [f"{insertion_files_path}.insertions.csv"]
    else:
        raise ValueError("Wrong path for insertion files")

    insname = ins_files[0].split("/")[-1].replace(".fasta.insertions.csv", "")

    ins_dict = {}  # {header: csv file} (only for sequences with good insertions)
    with open(f"{path_out}CL_squared_step_03/{insname}_sequences_with_bad_insertions.txt", "w") as doc:
        for ifile in ins_files:

            try:
                df = pd.read_csv(ifile, sep=None, engine="python", index_col="strain").fillna("")

                to_cons_cols = []   # ORFs
                bad_cols = []   # ORFs in START/STOP codon
                for col in df.columns:

                    # The genomic positions that are annotated in the csv file corresponds to the genomic positions of
                    # the reference sequence. If an insertion occurs before the genomic position 1, it is annotated as
                    # "@ref pos 0"
                    pos = int(col.split(" ")[-1])

                    for portion in portions.keys():
                        pos_list = ast.literal_eval(portions[portion])

                        # If I'm inspecting a genomic position comprised in a coding region (--> the output file
                        # annotates the genomic positions)
                        if pos in range(int(pos_list[0]), int(pos_list[1])):

                            # If I'm inspecting a genomic position that is not composed in a START or in a STOP codon
                            if pos not in range(int(pos_list[0]), int(pos_list[0])+2) and pos not in range(int(pos_list[1])-2, int(pos_list[1])):
                                to_cons_cols.append(col)
                                break

                            # If the START or the STOP codon is "broken"
                            else:
                                to_cons_cols.append(col)
                                bad_cols.append(col)

                for idx in df.index:
                    flag = 0

                    # Only the columns annotated before are inspected (genomic positions located into coding regions).
                    for col in to_cons_cols:

                        # The START or the STOP codon are "broken"
                        if col in bad_cols and len(df.at[idx, col]) != 0:
                            flag = 1
                            doc.write("{}\t{}\n".format(idx, ifile.split("\\")[-1]))
                            break

                        # Insertion in a coding region that is not divisible by 3
                        elif len(df.at[idx, col]) != 0 and len(df.at[idx, col]) % 3 != 0:
                            flag = 1
                            doc.write("{}\t{}\n".format(idx, ifile.split("\\")[-1]))
                            break

                    if flag == 0:
                        ins_dict[idx] = {}
                        ins_dict[idx]["file"] = ifile

                ins_df = pd.DataFrame.from_dict(ins_dict).T
                ins_df.to_csv(f"{path_out}CL_squared_step_03/{insname}_insertions_files.tsv", sep="\t")

            except KeyError:
                continue


def sequences_work(camp, meta, trim, ref_path, head, tail, path_out):

    all_sequences = list(meta.get_column("Headers"))
    print(f"ALL SEQUENCES: {all_sequences}")

    # ##############################################################################################################
    # Insertion files
    insname = camp.split("/")[-1].replace(".fasta", "")

    with open(f"{path_out}CL_squared_step_03/Sequence_cleaning_summary.txt", "a") as summary:

        summary.write(f"Working on file: {insname}\n")

        # Sequences with insertions that are non-divisible by 3 and are located in coding regions or with insertions in
        # the START or in the STOP codon
        bad_ins = []
        with open(f"{path_out}CL_squared_step_03/{insname}_sequences_with_bad_insertions.txt", "r") as doc:
            for line in doc:
                bad_ins.append(line.split("\t")[0])
        print(f"\n\nBAD INS: {bad_ins}")

        print(f"STARTING CLEANED SEQUENCES: {len(all_sequences)}")
        print(f"SEQUENCES WITH OUT OF FRAME INSERTIONS: {len(bad_ins)}")
        cleaned_seqs = list(set(all_sequences)-set(bad_ins))

        # Parameters to count ----------
        written_seqs = 0
        # ------------------------------

        alg = SeqIO.index(camp, "fasta")

        with open(f"{path_out}CL_squared_step_03/Cleaned_sequences_step_03/{insname}_cleaned_step_03.fasta", "w") as cl_fasta:
            print("ENTERED IN FASTA ----------")

            # Reference sequence
            if ref_path:
                reference = SeqIO.read(ref_path, "fasta")
                SeqIO.write(reference, cl_fasta, "fasta")

            for csec in cleaned_seqs:
                try:
                    sec = alg[csec]

                    #  INPUTS CONTROLLER ----------
                    if head:
                        head_trim = head
                    else:
                        head_trim = 1

                    if tail:
                        tail_trim = tail
                    # Tail assignment --> after seqs reading (length of each sequence --> no trimming)

                    if trim == "yes":

                        # Tail trimming
                        if tail:
                            sec.seq = sec.seq[0:tail_trim]
                        else:
                            sec.seq = sec.seq

                        # Head trimming
                        if head:
                            sec.seq = sec.seq[head_trim-1:]
                        else:
                            sec.seq = sec.seq

                    SeqIO.write(sec, cl_fasta, "fasta")
                    written_seqs += 1

                except KeyError:
                    print(f"Header {csec} not found in FASTA file")

        summary.write("STARTING LENGTH METADATA: {}\n".format(len(meta)))
        summary.write("TOT WRITTEN SEQS: {}\n".format(written_seqs))

        # Removal of metadata corresponding to discarded seqs ----------
        meta_new = meta.filter(pl.col("Headers").is_in(bad_ins))
        print("Final metadata length: {}\n".format(len(meta_new)))
        print(meta_new)

        cname = camp.split("/")[-1].replace(".fasta", "")
        print(f"{path_out}CL_squared_step_03/Temp_META_3/{cname}_step_03_BAD_INS.tsv")

        meta_new.write_csv(f"{path_out}CL_squared_step_03/Temp_META_3/{cname}_step_03_BAD_INS.tsv", separator="\t")

    return


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="The sequences are cleaned on the basis of the insertions that were "
                                             "trimmed by the augur align pipeline.\nSequences with out of frame "
                                             "insertions in coding regions are deleted.\nThe user can decide to only"
                                             " filter his sequences or, additionally, to also trim their extremities.",
                                             formatter_class=RawTextHelpFormatter)

# Pathogen coding regions
parser.add_argument("pathogen_portions_file", type=str, help="Path of the .tsv file containing the name of the "
                                                             "pathogen genome coding regions, each one associated "
                                                             "with the corresponding genomic positions. The file MUST "
                                                             "contain two columns: the first one named \"Name\" "
                                                             "containing the name of each coding region, the second one"
                                                             " named \"Positions\" containing the genomic positions of"
                                                             " each coding region.")

# augur align results folder (.fasta files and insertions .csv files)
parser.add_argument("alignment_files_path", type=str, help="Path of the folder containing the FASTA files "
                                                           "with the aligned sequences AND THE INSERTIONS CSV FILE. "
                                                           "ATTENTION: MAKE SURE THE NAME OF THE ALIGNMENT "
                                                           "FASTA FILE IS THE SAME AS THE NAME OF THE INSERTION .CSV "
                                                           "FILE WITH THE EXCEPTION OF THE FILE EXTENSION (.FASTA vs "
                                                           ".FASTA.INSERTIONS.CSV)")

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

# User metadata file
parser.add_argument("usr_metadata", type=str, help="Cleaned metadata file path.\nIf you have "
                                                   "multiple metadata files, concatenate or merge them.")

# Reference sequence path
parser.add_argument("--ref_path", type=str, help="Reference sequence path (optional).\nUse this parameter if you "
                                                 "want to insert the reference sequence as the first one of the "
                                                 "output FASTA files.")

# Sequence trimming
parser.add_argument("trim", type=str, help="Sequence trimming. If you want to trim the head and/or the tail of your "
                                           "sequences set this parameter as \"yes\", otherwise set it as \"no\". If you"
                                           " set it as yes pay attention that the default parameters of the head and "
                                           "the tail are head=0 and tail=length of your alignment sequences,so no "
                                           "trimming in both the head and the tail of the sequence.")

# Head trimming starting point
parser.add_argument("--head", type=int, help="Starting genomic position for the trimming of the head of the sequence"
                                             "(optional).\nDefault value: 0.")

# Tail trimming starting point
parser.add_argument("--tail", type=int, help="Starting genomic position for the trimming of the tail of the sequence."
                                             "(optional).\nDefault value: length of your sequences.")

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

# MAIN ####################
if __name__ == "__main__":

    pl.Config.set_tbl_formatting("ASCII_FULL_CONDENSED")

    if args.alignment_files_path.endswith(".fasta"):
        fastas = [args.alignment_files_path]
    elif args.alignment_files_path.endswith("/"):
        fastas = [i for i in glob.glob("{}*.fasta".format(args.alignment_files_path))]
    else:
        fastas = [i for i in glob.glob("{}/*.fasta".format(args.alignment_files_path))]

    print(fastas)
    print(len(fastas))

    # User metadata
    meta = pl.read_csv(args.usr_metadata, separator="\t", quote_char=None, ignore_errors=True).fill_null(value="Undefined")

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
    # GISAID ---------
    if type == "GISAID":

        with open(fastas[0], "r") as temp:
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
                                    + pl.col("Collection date") + "|"
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

            if not args.index_col_name:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the index column")

            if not args.header_format:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the header format")

            header_elements = args.header_format.split(",")

            inspected = 0

            all_idxs = list(meta.get_column(args.index_col_name))

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
            print("Headers ready")

    outfolder = args.outfolder
    if outfolder:
        if outfolder.endswith("/"):
            path_out = outfolder
        else:
            path_out = f"{outfolder}/"
    else:
        path_out = ""
    print(f"############## {path_out}")

    if not os.path.isdir(f"{path_out}CL_squared_step_03/"):
        os.makedirs(f"{path_out}CL_squared_step_03/")

    if not os.path.isdir(f"{path_out}CL_squared_step_03/Temp_META_3/"):
        os.makedirs(f"{path_out}CL_squared_step_03/Temp_META_3/")

    if not os.path.isdir(f"{path_out}CL_squared_step_03/Cleaned_sequences_step_03/"):
        os.makedirs(f"{path_out}CL_squared_step_03/Cleaned_sequences_step_03/")

    portions = portions_dict(args.pathogen_portions_file)

    # Insertion files ----------
    for camp in fastas:
        insertions_inspector(camp, portions, path_out)

    # MULTIPROCESSING -------------
    if os.path.isfile(f"{path_out}CL_squared_step_03/Sequence_cleaning_summary.txt"):
        os.remove(f"{path_out}CL_squared_step_03/Sequence_cleaning_summary.txt")

    multiprocessing.set_start_method('spawn')

    if len(fastas) == 1:
        sequences_work(fastas[0], meta, args.trim, args.ref_path, args.head, args.tail, path_out)

    else:
        processes = []
        for camp in fastas:

            p = multiprocessing.Process(target=sequences_work, args=(camp,
                                                                    meta,
                                                                    args.trim,
                                                                    args.ref_path,
                                                                    args.head,
                                                                    args.tail,
                                                                    path_out))
            processes.append(p)

        for p in processes:
            p.start()

        for process in processes:
            process.join()

    all_mini_meta = [i for i in glob.glob(f"{path_out}CL_squared_step_03/Temp_META_3/*.tsv")]

    new_mm = pd.read_csv(all_mini_meta[0], sep="\t")
    print(new_mm)
    all_meta_seqs = len(new_mm)
    for mm in all_mini_meta[1:]:
        meta_file = pd.read_csv(mm, sep="\t")
        all_meta_seqs += len(meta_file)
        new_mm = pd.concat([new_mm, meta_file], sort=False)

    bad_head = list(new_mm["Headers"])

    cleaned_meta = meta.filter(~pl.col("Headers").is_in(bad_head))
    print(len(cleaned_meta))

    cleaned_meta.write_csv(f"{path_out}CL_squared_step_03/Cleaned_meta_step_03.tsv", separator="\t")

    shutil.rmtree(f"{path_out}CL_squared_step_03/Temp_META_3")
