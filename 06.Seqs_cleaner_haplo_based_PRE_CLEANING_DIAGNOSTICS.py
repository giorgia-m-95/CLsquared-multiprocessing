
# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - User FASTA files to be cleaned:
#
# - Out folder for the summary file and the diagnostics .tsv files (optional).
#   Default out folder: current folder.
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# The script inspects the sequences that the user wants to filter in order to provide useful information to consciously
# clean the input dataset. Information related to the parameters that can be chosen in the cleaning step 06 will be
# contained in the output .tsv file.
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - Summary .txt file
# - Diagnostics .tsv files
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

# FUNCTIONS ############################################################################################################


def sequences_inspection(camp, tsvs, usr_metadata, path_out):

    fasta_name = camp.split("/")[-1].replace(".fasta", "")
    print(f"Working on {fasta_name} ----------")

    meta = pl.read_csv(usr_metadata, separator='\t', quote_char=None, ignore_errors=True).fill_null(value="Undefined")

    pl.Config.set_tbl_formatting("ASCII_FULL_CONDENSED")
    print(meta)

    with open(f"{path_out}CL_squared_PRE_step_06-DIAGNOSTICS/{fasta_name}_sequences_inspection_diagnostics_summary.txt", "w") as summary:

        for tsv in tsvs:
            if fasta_name in tsv:
                print("ENTERED")

                df = pl.read_csv(tsv, separator='\t', ignore_errors=True).fill_null(value="Undefined")

                # We work on the defining mutations (clustering croup 1)
                df = df.filter(pl.col("Group") == 1)

                # Defining haplotype definition for the clade we are working on
                haplo = {}
                indexes = sorted(list(df.get_column("Genomic_position")))
                for idx in indexes:
                    haplo[idx] = {}
                    value = df.filter(pl.col("Genomic_position") == idx).select("ref").item()
                    haplo[idx]["ref"] = value

                    major_el = 0
                    major_el_col = ""
                    for col in list(df.columns)[2:6]:
                        value = df.filter(pl.col("Genomic_position") == idx).select(col).item()
                        value = float(value)
                        if value > major_el:
                            major_el = value
                            major_el_col = col

                    haplo[idx]["mut"] = major_el_col

                # Sequences inspection (the possibility to have multiple FASTA files for each clade is taken into account)
                ffiles = [camp]

                all_headers = list(meta.get_column("Headers"))

                seqs_dict = {}
                for ffile in ffiles:
                    name = ffile.split("/")[-1].replace(".fasta", "")
                    tot_seqs = 0
                    summary.write(f"Inspecting file {ffile}\n")

                    seqs = SeqIO.index(ffile, "fasta")

                    summary.write(f"TOT SEQUENCES TO INSPECT: {len(seqs)}\n")

                    for head in all_headers:
                        try:
                            sec = seqs[head]
                            tot_seqs += 1
                            wild = 0
                            diff = 0
                            for pos in haplo.keys():
                                idx = pos - 1
                                if sec.seq[idx] == haplo[pos]["ref"]:
                                    wild += 1
                                elif sec.seq[idx] != haplo[pos]["mut"]:
                                    diff += 1
                            seqs_dict[sec.id] = {}
                            seqs_dict[sec.id]["Reversions"] = wild
                            seqs_dict[sec.id]["Other mutations"] = diff
                            seqs_dict[sec.id]["Diff percentage"] = f"{(wild+diff)/len(haplo.keys())*100}%"
                        except KeyError:
                            continue

                    summary.write(f"TOT INSPECTED SEQUENCES: {len(seqs)}\n\n")

                    df_to_save = pd.DataFrame.from_dict(seqs_dict).T
                    df_to_save.to_csv(f"{path_out}CL_squared_PRE_step_06-DIAGNOSTICS/"
                                      f"{name}_DIAGNOSTICS.tsv", sep="\t")

    return


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="The script inspects the sequences that the user wants to filter in order "
                                             "to provide useful information to consciously clean the input dataset. "
                                             "Information related to the parameters that can be chosen in the cleaning "
                                             "step 06 will be contained in the output .tsv file.",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("usr_fastas", type=str, help="Path of the FASTA file that the user wants to filter in "
                                                 "the next CLsquared step (06).\nIf you want to "
                                                 "inspect multiple FASTA files, please insert the folder that"
                                                 "contains all your files.\nATTENTION: IF THE CHOSEN FOLDER"
                                                 " CONTAINS EXTRA FASTA FILES, THEY WILL BE INSPECTED BY THE SCRIPT."
                                                 " PLEASE VERIFY THAT THE NAME OF THE FASTA FILE IS CONTAINED IN THE "
                                                 "NAME OF THE CLUSTERING RESULT FILE.")

# Clustering step results path
parser.add_argument("clustering_res_path", type=str, help="Path of the .tsv files containing the results of the "
                                                          "clustering step (step 05 of CLsquared pipeline, "
                                                          "Clustered_prevalences folder).\nATTENTION: THE NAME OF THE "
                                                          "ALIGNED FASTA FILE TO BE CLEANED SHOULD BE CONTAINED IN THE "
                                                          "NAME OF THE CLUSTERING RESULTS FILE.")

# User metadata file
parser.add_argument("usr_metadata", type=str, help="User metadata file's path.\nThe metadata have to be correctly "
                                                   "associated to the sequences that you want to inspect.\nIf you have "
                                                   "multiple metadata files, concatenate or merge them.")
# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

args = parser.parse_args()

if __name__ == "__main__":

    # INPUTS CONTROLLER ----------
    if args.usr_fastas.endswith(".fasta"):
        camps = [args.usr_fastas]
    elif args.usr_fastas.endswith("/"):
        camps = [i for i in glob.glob(f"{args.usr_fastas}*.fasta")]
    else:
        camps = [i for i in glob.glob(f"{args.usr_fastas}/*.fasta")]

    if args.usr_metadata.endswith(".tsv"):
        all_meta = [args.usr_metadata]
    elif args.usr_fastas.endswith("/"):
        all_meta = [i for i in glob.glob(f"{args.usr_metadata}*.tsv")]
    else:
        all_meta = [i for i in glob.glob(f"{args.usr_metadata}/*.tsv")]

    if args.clustering_res_path.endswith(".tsv"):
        tsvs = [args.clustering_res_path]
    elif args.clustering_res_path.endswith("/"):
        tsvs = [i for i in glob.glob(f"{args.clustering_res_path}*.tsv")]
    else:
        tsvs = [i for i in glob.glob(f"{args.clustering_res_path}/*.tsv")]

    if args.outfolder:
        if args.outfolder.endswith("/"):
            path_out = args.outfolder
        else:
            path_out = f"{args.outfolder}/"
    else:
        path_out = ""

    if not os.path.isdir(f"{path_out}CL_squared_PRE_step_06-DIAGNOSTICS"):
        os.makedirs(f"{path_out}CL_squared_PRE_step_06-DIAGNOSTICS")

    multiprocessing.set_start_method('spawn')

    processes = []
    if len(camps) == 1:
        sequences_inspection(camps[0], tsvs, args.usr_metadata, path_out)

    else:
        for camp in camps:

            p = multiprocessing.Process(target=sequences_inspection, args=(camp,
                                                                           tsvs,
                                                                           args.usr_metadata,
                                                                           path_out))

            processes.append(p)

        for p in processes:
            p.start()

        for process in processes:
            process.join()
