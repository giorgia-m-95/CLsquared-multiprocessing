# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - User FASTA files, aligned VS the pathogen reference sequence, to be cleaned
#
# - User's metadata.
#   MANDATORY MINIMUM FORMAT:
#    - columns: Accession ID (sequence ID) (GISAID metadata format)
#               Accession (sequence ID) (NCBI metadata format)
#               If your metadata file doesn't match the GISAID and NCBI formats you can provide the name of the column
#		containing the SEQUENCE ID THAT IS THE ONLY MANDATORY ONE. If it is unnamed, use its index number
#		(e.g., 0, 1, ...).
#		You can also provide a column containing the headers of the sequences: "Headers" column (see previuos
#               scripts of the pipeline)
#
#    - metadata have to be associated to the sequences contained in the FASTA files
#
# - Out folder for the summary file and the cleaned metadata and FASTA files (optional).
#   Default out folder: current folder.
# ----------------------------------------------------------------------------------------------------------------------

# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# - FASTA aligned sequence filtering.
#   If required by the user, the sequences will be discarded on the basis of the following criteria:
#    - A percentage of ambiguous characters
#    - Head/tail gap trails longer than x bases (default: 90) # TODO MANTENERE QUESTA SOGLIA DI MINIMO?
#    - Gaps causing frameshifting
# - Metadata cleaning (metadata linked to deleted sequences are removed)
# ----------------------------------------------------------------------------------------------------------------------

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
from Bio import AlignIO, SeqIO
import glob
import os
import pandas as pd
import polars as pl
import re


# FUNCTIONS ############################################################################################################
# Pathogen coding regions ----------
def portions_dict(pathogen_portions_file):

    # PANDAS ----------
    # portions_df = pd.read_csv(pathogen_portions_file, sep=None, engine="python")
    # --------------------

    portions_df = pl.read_csv(pathogen_portions_file, separator='\t', ignore_errors=True).fill_null(value="Undefined")

    portions = {}

    # PANDAS ----------
    # for idx in portions_df.index:
    #     portions[portions_df.at[idx, "Name"]] = ast.literal_eval(portions_df.at[idx, "Positions"])
    # --------------------

    portion_names = portions_df.get_column("Name")

    if len(portion_names) != len(list(set(portion_names))):
        raise ValueError("ATTENTION! Duplicate portion names")
    else:
        for pn in portion_names:
            portions[pn] = portions_df.select(pl.col("Name") == pn, pl.col("Positions"))

    return portions


def aligned_seqs_cleaner(camp, usr_meta, ambiguous, tolerance, trails, head, tail, frameshifting, portions, outfolder, index_col_name):

    # Metadata controller ----------
    meta = pl.read_csv(usr_meta, separator="\t", ignore_errors=True).fill_null(value="Undefined")

    # Input columns checking ##########
    meta_columns = meta.columns
    # Public metadata
    if not index_col_name:
        if "Accession ID" not in meta_columns and "Accession" not in meta_columns:
            raise ValueError(
                "Wrong column names or missing columns. Minimum mandatory column for GISAID or NCBI metadata:"
                " Sequence accession ID")
        # GISAID metadata
        elif "Accession ID" in meta_columns:
            type = "GISAID"
        # NCBI metadata
        elif "Accession" in list(meta.columns):
            meta.set_index("Accession", inplace=True)
            type = "NCBI"
    # Personal metadata
    else:
        if index_col_name.isdigit():
            meta.set_index(int(index_col_name), inplace=True)
        else:
            meta.set_index(index_col_name, inplace=True)
        type = "personal"

    with open("{}CL_squared_step_04/Sequence_cleaning_summary.txt".format(path_out), "w") as summary:

	# Parameters to count ----------
        disc_ambiguous = 0
        disc_head_tail = 0
        disc_out_frame = 0
        epi_disc = []
        file_counter = 0
        in_seqs = 0
        tot_written_seqs = 0
        # ---------------------

        # Headers making ##########
        # GISAID ---------
        if type == "GISAID":
            print(file)
            with open(file, "r") as temp:
                for line in temp:
                    if "EPI" in line:
                        type = "A"
                        print(type)
                    else:
                        type = "B"
                        print(type)
                    break

            if type == "B" and "Headers" not in meta.columns:
                meta = meta.with_columns(pl.concat_str([pl.col("Virus name"), pl.col("Collection date"),
                                         pl.col("Submission date")], separator="|").alias("Headers"))
                meta.with_columns(pl.col("Headers").str.replace("|Headers", ""))

            elif type == "A" and "Headers" not in meta.columns:
                meta = meta.with_columns(pl.concat_str([pl.col("Virus name"), pl.col("Accession ID"),
                                         pl.col("Collection date")], separator="|").alias("Headers"))
                meta.with_columns(pl.col("Headers").str.replace("|Headers", ""))

            print("Headers ready")

        # NCBI ----------
        elif type == "NCBI":
            if "Headers" not in meta.columns:
                meta = meta.with_columns(pl.get_column(index_col).alias("Headers"))
                print("Headers ready")

        # Personal ----------
        elif type == "personal":
            if "Headers" not in meta.columns:
                meta = meta.with_columns(pl.get_column(index_col).alias("Headers"))
                print("Headers ready")

        # SEQUENCES PARSING AND CLEANING ##########
        headers = meta.get_column("Headers")

        all_msa = [camp] # -------------------> attenzione, ricontrollare
        for msa in all_msa:
            file_counter += 1

            algn = AlignIO.index(msa, "fasta")
            summary.write("{}: {} starting sequences\n".format(msa.split("/")[-1], len(algn)))

            with open("{}CL_squared_step_04/FASTA_file_{}_cleaned_step_04.fasta".format(path_out, file_counter), "w") as doc:

                for header in headers:

		    # INDENTARE ------------------------------------------------------------
		    try:
                        sec = algn[header].seq
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
                                        epi_disc.append(record.id.split("/")[1])
                                    else:
                                        # If required, gaps causing frameshift check
                                        if frameshifting == "yes":
                                            out_flag = frameshifting_checker(sec, portions)
                                            if out_flag == 1:
                                                epi_disc.append(record.id.split("/")[1])
                                                disc_out_frame += 1
                                            elif out_flag == 0:
                                                SeqIO.write(record, doc, "fasta")
                                                tot_written_seqs += 1
                                        else:
                                            SeqIO.write(record, doc, "fasta")
                                            tot_written_seqs += 1
                                else:
                                    if frameshifting == "yes":
                                        out_flag = frameshifting_checker(sec, portions)
                                        if out_flag == 1:
                                            epi_disc.append(record.id.split("/")[1])
                                            disc_out_frame += 1
                                        elif out_flag == 0:
                                            SeqIO.write(record, doc, "fasta")
                                            tot_written_seqs += 1
                                    else:
                                        SeqIO.write(record, doc, "fasta")
                                        tot_written_seqs += 1
                            else:
                                epi_disc.append(record.id.split("/")[1])
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
                                    epi_disc.append(record.id.split("/")[1])
                                else:
                                    if frameshifting == "yes":
                                        out_flag = frameshifting_checker(sec, portions)
                                        if out_flag == 1:
                                            epi_disc.append(record.id.split("/")[1])
                                            disc_out_frame += 1
                                        elif out_flag == 0:
                                            SeqIO.write(record, doc, "fasta")
                                            tot_written_seqs += 1
                                    else:
                                        SeqIO.write(record, doc, "fasta")
                                        tot_written_seqs += 1
                            else:
                                if frameshifting == "yes":
                                    out_flag = frameshifting_checker(sec, portions)
                                    if out_flag == 1:
                                        epi_disc.append(record.id.split("/")[1])
                                        disc_out_frame += 1
                                    elif out_flag == 0:
                                        SeqIO.write(record, doc, "fasta")
                                        tot_written_seqs += 1
                                else:
                                    SeqIO.write(record, doc, "fasta")
                                    tot_written_seqs += 1

        # Print reports ----------
        summary.write("\nTOT STARTING SEQUENCES: {}\n".format(in_seqs))
        summary.write("TOT WRITTEN SEQUENCES: {}\n".format(tot_written_seqs))
        summary.write("TOT DISCARDED SEQUENCES\n: {}".format(len(epi_disc)))
        if trails == "yes":
            summary.write('SEQS WITH HEAD-TAIL GAPS: {}\n'.format(disc_head_tail))
        if frameshifting == "yes":
            summary.write('FRAMESHIFT SEQS: {}\n'.format(disc_out_frame))

        # Removal of metadata corresponding to discarded seqs ----------

        # PANDAS ----------
        # meta = pd.read_csv(usr_meta, sep=None, engine="python", index_col="Accession ID")
        # --------------------

        # summary.write("\nStarting metadata length: {}\n".format(len(meta)))

        # PANDAS ----------
        # meta = meta[~meta.index.isin(epi_disc)]
        # --------------------

        meta = meta.filter(~pl.col("Hreaders").is_in(epi_disc))
        print("Final metadata length: {}\n".format(len(meta)))

        if len(meta) != tot_written_seqs:
            raise ValueError("ERROR! Final metadata do not correspond to written sequences")

	cname = camp.replace(".fasta", "")
        meta.to_csv("{}CL_squared_step_04/Temp_META_4/{cname}_step_04.tsv".format(path_out), sep="\t")


# Sequences containing gaps causing frameshift are discarded ----------
def frameshifting_checker(sec, portions):
    flag = 0

    for orf in portions:
        start = portions[orf][0] - 1
        stop = portions[orf][1]
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
                                             "default: 90)\n- Gaps causing frameshift in the coding regions\nMetadata "
                                             "cleaning: metadata linked to deleted sequences are removed",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="User FASTA file's path.\nIf you want to inspect multiple FASTA files,"
                                                 " write the name of the folder that contains all your files.\nEvery "
                                                 "sequence has to be associated to the corresponding metadata "
                                                 "(usr_meta).\nATTENTION: IF THE CHOSEN FOLDER CONTAINS EXTRA FASTA "
                                                 "FILES, THEY WILL BE INSPECTED BY THE SCRIPT.\nMandatory sequence's "
                                                 "header format: >Place/Accession ID/Collection date.\nIf the format of"
                                                 " your headers does not match this one, a formatting script is "
                                                 "provided:\nHeader_formatter.py")

# User metadata file
parser.add_argument("usr_meta", type=str, help="User metadata file's path.\nThe metadata have to be correctly "
                                               "associated to the sequences that you want to inspect.\nIf you have "
                                               "multiple metadata files, concatenate or merge them.")

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
                                                    "frameshifts in the coding regions.\n")

# Pathogen coding regions
parser.add_argument("pathogen_portions_file", type=str, help="Path of the .csv/.tsv file containing the name of the "
                                                             "pathogen genome coding regions, each one associated "
                                                             "with the corresponding genomic positions. The file MUST "
                                                             "contain two columns: the first one named \"Name\" "
                                                             "containing the name of each coding region, the second one"
                                                             " named \"Positions\" containing the genomic positions of"
                                                             " each coding region.")

# Index column name (sequence id)
parser.add_argument("--index_col_name", type=str, help="If your metadata file doesn't match the GISAID and NCBI "
                                                        "formats you can provide the name of the column containing the "
                                                        "SEQUENCE ID THAT IS THE ONLY MANDATORY ONE. If it is unnamed, "
                                                        "use its index number (e.g., 0, 1, ...).")

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")


#  INPUTS CONTROLLER ----------
if outfolder:
    if outfolder.endswith("/"):
        path_out = outfolder
    else:
        path_out = "{}/".format(outfolder)
else:
    path_out = ""

if not os.path.isidr(f"{path_out}CL_squared_step_04/Temp_META_4/"):
    os.makedirs(f"{path_out}CL_squared_step_04/Temp_META_4/")

if frameshifting not in ["yes", "no"]:
    raise ValueError("You have entered the wrong input. User-allowed inputs: yes, no")

if ambiguous not in ["yes", "no"]:
    raise ValueError("You have entered the wrong input. User-allowed inputs: yes, no")

if ambiguous == "no" and tolerance:
    raise ValueError("You have set \"ambiguous\" parameter as \"no\": set ambiguous parameter as \"yes\" if you "
                     "want to filter your sequences on the basis of the ambiguous characters.")

if trails == "no" and head or trails == "no" and tail:
    raise ValueError("You have set \"trail\" parameter as \"no\": set \"trail\" parameter as \"yes\" if you want to"
                         " filter your sequences on the basis of the head/tail gap trails.")

if usr_fastas.endswith(".fasta"):
    all_msa = [usr_fastas]
elif usr_fastas.endswith("/"):
    all_msa = [i for i in glob.glob("{}*.fasta".format(usr_fastas))]
else:
    all_msa = [i for i in glob.glob("{}/*.fasta".format(usr_fastas))]

if not os.path.isdir("{}CL_squared_step_04".format(path_out)):
    os.makedirs("{}CL_squared_step_04".format(path_out))


# Outputs path ----------

args = parser.parse_args()

portions = portions_dict(args.pathogen_portions_file)

aligned_seqs_cleaner(args.usr_fastas, args.usr_meta, args.ambiguous, args.tolerance, args.trails, args.head, args.tail,
                     args.frameshifting, portions, args.outfolder)


# MULTIPROCESSING -------------

processes = []
for camp in all_msa:
    p = multiprocessing.Process(target=sequence_cleaner, args=(camp,
                                                                   args.usr_meta,
                                                                   args.ambiguous,
                                                                   args.tolerance,
                                                                   args.trails,
                                                                   args.head,
                                                                   args.tail,
                                                                   args.frameshifting,
                                                                   args.portions,
                                                                   path_out,
                                                                   args.index_col_name))
    processes.append(p)
    p.start()

for process in processes:
    process.join()

all_mini_meta = [i for i in glob.glob(f"{path_out}CL_squared_step_04/Temp_META_2/*.tsv")]

new_mm = pd.read_csv(all_mini_meta[0], sep="\t")
print(new_mm)
all_meta_seqs = len(new_mm)
for mm in all_mini_meta[1:]:
    meta_file = pd.read_csv(mm, sep="\t")
    all_meta_seqs += len(meta_file)
    new_mm = pd.concat([new_mm, meta_file], sort=False)
print(len(new_mm))

new_mm.to_csv(f"{path_out}CL_squared_step_04/Cleaned_meta_step_02_2.tsv", sep="\t")

shutil.rmtree(f"{path_out}CL_squared_step_04/Temp_META_2")


