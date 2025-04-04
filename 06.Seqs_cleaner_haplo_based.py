# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - .tsv files containing the results of the clustering step
#
# - FASTA files containing the msa of your sequences aligned versus the pathogen reference sequence
#
# - FASTA files containing the "not-aligned" sequences if you want to work with them
#
# - User metadata
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# This script is designed to clean your sequences on the basis of the clustering step 05.
# The sequences are evaluated on the basis of the comparison between their pattern of mutations and the defining
# haplotype computed in the clustering step.
# The clustering step is meant to be an additional cleaning step and not a new way to call the defining mutations of a
# specific viral variant.
# The user can choose the stringency of the cleaning step through a tolerance threshold for the percentage of admitted
# reversions (or different mutations) with respect to the novel defining haplotype.
# Reversions are the most evaluated mutations because of their inflated presence in GISAID due to sequencing/assembling
# errors.
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
from Bio import SeqIO
import glob
import multiprocessing
import os
import pandas as pd
import polars as pl


# FUNCTIONS ############################################################################################################

def seqs_cleaner(tsvs, camp, usr_metadata, threshold, path_out):

    name = camp.split("/")[-1].replace(".fasta", "")
    print(f"Working on file {name} ----------")

    with open(f"{path_out}CL_squared_step_06/{name}_msa_and_meta_cleaning_summary.txt", "w") as summary:
        for tsv in tsvs:

            if name in tsv:

                meta = pl.read_csv(usr_metadata, separator='\t', quote_char=None, ignore_errors=True).fill_null(value="Undefined")

                # Recombinant clade is not taken into account in our analyses
                summary.write(f"Entered file {name}\n")

                df = pl.read_csv(tsv, separator='\t', ignore_errors=True).fill_null(value="Undefined")

                # We work on the defining mutations (clustering croup 1)
                df = df.filter(pl.col("Group") == 1)

                # Defining haplotype definition for the clade we are working on
                haplo = {}
                indexes = list(df.get_column("Genomic_position"))
                indexes = [int(i) for i in indexes]
                for idx in sorted(indexes):
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

                haplo_df = pd.DataFrame.from_dict(haplo).T
                haplo_df.to_csv(
                    f"{path_out}CL_squared_step_06/Defining_haplotypes/{name}_defining_haplo.tsv",
                    sep="\t")

                ffiles = [camp]

                for ffile in ffiles:
                    seqs = SeqIO.index(ffile, "fasta")

                    summary.write(f"Tot starting seqs {name}: {len(seqs)}\n")

                    all_headers = list(meta.get_column("Headers"))

                    with open(
                            f"{path_out}CL_squared_step_06/Cleaned_msa_seqs/{name}_haplo_cleaned.fasta",
                            "w") as cleaned_fasta:
                        with open(f"{path_out}CL_squared_step_06/Discarded_msa_seqs/{name}_haplo_cleaned_disc.fasta",
                                  "w") as discarded_fasta:

                            disc_epis = []
                            surviving = 0
                            discarded = 0
                            for head in all_headers:

                                try:
                                    sec = seqs[head]

                                    diff = 0
                                    for pos in haplo.keys():
                                        idx = pos - 1
                                        if sec.seq[idx] == haplo[pos]["mut"]:
                                            continue
                                        else:
                                            diff += 1

                                    if threshold:
                                        if (diff / len(haplo.keys()))*100 < threshold:
                                            sec.description = ""
                                            SeqIO.write(sec, cleaned_fasta, "fasta")
                                            surviving += 1
                                        else:
                                            disc_epis.append(sec.id)
                                            SeqIO.write(sec, discarded_fasta, "fasta")
                                            discarded += 1
                                    else:
                                        if diff == 0:
                                            sec.description = ""
                                            SeqIO.write(sec, cleaned_fasta, "fasta")
                                            surviving += 1
                                        else:
                                            disc_epis.append(sec.id)
                                            SeqIO.write(sec, discarded_fasta, "fasta")
                                            discarded += 1

                                except KeyError:
                                    continue

                print(f"Discarded sequences: {discarded}")
                print(f"Surviving sequences: {surviving}")

                if discarded + surviving != len(seqs):
                    raise ValueError("ERROR")

                # Metadata cleaning ----------
                start_length = len(meta)
                summary.write(f"STARTING LENGTH META: {start_length}\n")

                with open(f"{path_out}CL_squared_step_06/Discarded_msa_seqs/{name}_discarded_epis.txt", "w") as disc_file:
                    for de in disc_epis:
                        disc_file.write(f"{de}\n")

                summary.write(f"TOTAL SEQUENCES TO DISCARD: {len(disc_epis)}\n")

                meta_filtered = meta.filter(~pl.col("Headers").is_in(disc_epis))

                end_length = len(meta_filtered)
                summary.write(f"FINAL LENGTH META: {end_length}\n\n")

                if (len(disc_epis) + end_length) != start_length:
                    raise ValueError("ERROR!")

                meta_filtered.write_csv(f"{path_out}CL_squared_step_06/{name}_post_haplo_cleaning_meta.tsv", separator="\t")

                break

    return


# Trimmed sequences cleaning ----------
def trimmed_seqs_cleaner(trimmed_seq, meta, path_out):

    with open(f"{path_out}CL_squared_step_06/Deligned_sequences_cleaning_summary.txt", "w") as summary:

        discarded = 0
        written = 0
        starting_seqs = 0

        name = trimmed_seq.split("/")[-1].replace(".fasta", "")

        summary.write(f"Working on file {name}\n")

        seqs = SeqIO.index(trimmed_seq, "fasta")

        starting_seqs += len(seqs)

        with open(f"{path_out}CL_squared_step_06/Cleaned_dealigned_seqs/{name}_haplo_cleaned.fasta", "w") as cleaned_fasta:
            with open(f"{path_out}CL_squared_step_06/Discarded_dealigned_seqs/{name}_haplo_cleaned.fasta",
                      "w") as discarded_fasta:

                headers = list(meta.get_column("Headers"))

                for head in headers:
                    try:
                        sec = seqs[head]
                        sec.description = ""
                        SeqIO.write(sec, cleaned_fasta, "fasta")
                        written += 1

                    except KeyError:
                        sec.description = ""
                        SeqIO.write(sec, discarded_fasta, "fasta")
                        discarded += 1
                        continue

        summary.write(f"STARTING SEQUENCES NUMBER: {starting_seqs}\n")
        summary.write(f"FINAL SEQUENCES NUMBER: {written}\n")
        summary.write(f"DISCARDED SEQUENCES NUMBER: {discarded}")


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="This script cleans the sequences on the basis of the clustering step 05. "
                                             "The sequences are evaluated on the basis of the comparison between their "
                                             "pattern of mutations and the haplotype computed by the clustering step.",
                                 formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("msa_clades_path", type=str, help="FASTA files path.\nPlease insert the folder "
                                                      "containing all your files.\nATTENTION: IF THE FOLDER"
                                                      " CONTAINS EXTRA FASTA FILES, THEY WILL BE INSPECTED BY THE "
                                                      "SCRIPT.\nThe script cleans the sequences that are aligned versus"
                                                      " the reference sequence that has been used for the clustering "
                                                      "step.")

# "Dealigned" sequences FASTA files path
parser.add_argument("--trimmed_clades_path", type=str, help="Path of the FASTA files containing the "
                                                            "sequences that are not aligned versus the reference "
                                                            "sequence.\nPlease insert the folder containing all your "
                                                            "files.\nATTENTION: IF THE FOLDER CONTAINS EXTRA FASTA "
                                                            "FILES, THEY WILL BE INSPECTED BY THE SCRIPT.")

# Clustering step results path
parser.add_argument("clustering_res_path", type=str, help="Path of the .tsv file containing the results of "
                                                          "the clustering step 05 (Clustered prevalences folder). "
                                                          "ATTENTION: THE NAME OF THE FASTA FILE TO BE CLEANED SHOULD "
                                                          "BE CONTAINED IN THE NAME OF THE CLUSTERING RESULTS FILE.")

# User metadata file
parser.add_argument("usr_metadata", type=str, help="Metadata file path.\nIf you have multiple metadata "
                                                   "files, concatenate or merge them.")

# Tolerance threshold
parser.add_argument("--threshold", type=float, help="Tolerance threshold for the number of reversions or "
                                                    "different mutations with respect to the computed defining "
                                                    "haplotype (percentage value, optional). Default value: 0%%")

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: "
                                                  "current folder.")

args = parser.parse_args()

if __name__ == "__main__":

    #  INPUTS CONTROLLER ----------
    if args.clustering_res_path.endswith(".tsv"):
        tsvs = [args.clustering_res_path]
    elif args.clustering_res_path.endswith("/"):
        tsvs = [i for i in glob.glob(f"{args.clustering_res_path}*.tsv")]
    else:
        tsvs = [i for i in glob.glob(f"{args.clustering_res_path}/*.tsv")]

    if args.msa_clades_path.endswith("/"):
        msa_path = args.msa_clades_path
        all_msa = [i for i in glob.glob(f"{msa_path}*.fasta")]
    elif args.msa_clades_path.endswith(".fasta"):
        msa_path = args.msa_clades_path
        all_msa = [msa_path]
    else:
        msa_path = f"{args.msa_clades_path}/"
        all_msa = [i for i in glob.glob(f"{msa_path}*.fasta")]

    if args.outfolder:
        if args.outfolder.endswith("/"):
            path_out = args.outfolder
        else:
            path_out = f"{args.outfolder}/"
    else:
        path_out = ""

    if not os.path.isdir(f"{path_out}CL_squared_step_06"):
        os.makedirs(f"{path_out}CL_squared_step_06")

    # Cleaned sequences folder ----------
    if not os.path.isdir(f"{path_out}CL_squared_step_06/Cleaned_msa_seqs"):
        os.makedirs(f"{path_out}CL_squared_step_06/Cleaned_msa_seqs")

    # Discarded sequences folder ----------
    if not os.path.isdir(f"{path_out}CL_squared_step_06/Discarded_msa_seqs"):
        os.makedirs(f"{path_out}CL_squared_step_06/Discarded_msa_seqs")

    # Defining haplotype files
    if not os.path.isdir(f"{path_out}CL_squared_step_06/Defining_haplotypes"):
        os.makedirs(f"{path_out}CL_squared_step_06/Defining_haplotypes")


    multiprocessing.set_start_method('spawn')

    processes = []

    if len(all_msa) == 1:

        seqs_cleaner(tsvs, all_msa[0], args.usr_metadata, args.threshold, path_out)

        if args.trimmed_clades_path:

            if args.trimmed_clades_path.endswith(".fasta"):
                all_trimmed_seqs = [args.trimmed_clades_path]
            elif args.trimmed_clades_path.endswith("/"):
                all_trimmed_seqs = [i for i in glob.glob(f"{args.trimmed_clades_path}*.fasta")]
            else:
                all_trimmed_seqs = [i for i in glob.glob(f"/{args.trimmed_clades_path}*.fasta")]

            # Cleaned sequences folder ----------
            if not os.path.isdir(f"{path_out}CL_squared_step_06/Cleaned_dealigned_seqs"):
                os.makedirs(f"{path_out}CL_squared_step_06/Cleaned_dealigned_seqs")

            # Discarded sequences folder ----------
            if not os.path.isdir(f"{path_out}CL_squared_step_06/Discarded_dealigned_seqs"):
                os.makedirs(f"{path_out}CL_squared_step_06/Discarded_dealigned_seqs")

            cleaned_meta = f"{path_out}CL_squared_step_06/Post_haplo_cleaning_meta.tsv"

            trimmed_seqs_cleaner(all_trimmed_seqs, cleaned_meta, path_out)

    else:
        for camp in all_msa:

            p = multiprocessing.Process(target=seqs_cleaner, args=(tsvs,
                                                                   camp,
                                                                   args.usr_metadata,
                                                                   args.threshold,
                                                                   path_out))

            processes.append(p)

        for p in processes:
            p.start()

        for process in processes:
            process.join()

        if args.trimmed_clades_path:

            if args.trimmed_clades_path.endswith(".fasta"):
                all_trimmed_seqs = [args.trimmed_clades_path]
            elif args.trimmed_clades_path.endswith("/"):
                all_trimmed_seqs = [i for i in glob.glob(f"{args.trimmed_clades_path}*.fasta")]
            else:
                all_trimmed_seqs = [i for i in glob.glob(f"/{args.trimmed_clades_path}*.fasta")]

            # Cleaned sequences folder ----------
            if not os.path.isdir(f"{path_out}CL_squared_step_06/Cleaned_dealigned_seqs"):
                os.makedirs(f"{path_out}CL_squared_step_06/Cleaned_dealigned_seqs")

            # Discarded sequences folder ----------
            if not os.path.isdir(f"{path_out}CL_squared_step_06/Discarded_dealigned_seqs"):
                os.makedirs(f"{path_out}CL_squared_step_06/Discarded_dealigned_seqs")

            cleaned_meta = f"{path_out}CL_squared_step_06/Post_haplo_cleaning_meta.tsv"

            if len(all_trimmed_seqs) == 1:

                trimmed_seqs_cleaner(all_trimmed_seqs[0], cleaned_meta, path_out)

            else:
                for trimmed_seq in all_trimmed_seqs:

                    trimmed_seqs_cleaner(trimmed_seq, cleaned_meta, path_out)

                    p = multiprocessing.Process(target=trimmed_seqs_cleaner, args=(trimmed_seq,
                                                                                   cleaned_meta,
                                                                                   path_out))

                    processes.append(p)

                for p in processes:
                    p.start()

                for process in processes:
                    process.join()
