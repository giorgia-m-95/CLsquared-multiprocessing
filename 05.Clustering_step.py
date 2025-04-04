# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - User FASTA files, with sequences divided on the basis of a grouping criterion: e.g. the viral clade.
#   You can divide your sequences using the script C.1.
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
# The script has the objective to recognize the prevalent defining haplotype of a group of sequences starting from a
# cleaned dataset.
#
# WE SUGGEST YOU TO USE THIS SCRIPT AFTER EXPLOITING THE WHOLE PIPELINE ON YOUR DATA WITH THE MOST STRINGENT AS POSSIBLE
# FILTERS.
#
# You can personalize the clustering step as you prefer. The user can chose:
#   - the head/tail genomic positions to not be considered in the analysis
#   - A mutation-frequency threshold in order to select which frequencies to discard in the clustering analysis
#     (putative noise).
#
# The sequences can be further cleaned in the following steps of the pipeline on the basis of the clustering results.
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# Multiple .tsv and .html files reporting the results of the clustering of the mutation percentage for each clade.
# - Summary .txt file.
# ----------------------------------------------------------------------------------------------------------------------
#
# ######################################################################################################################


import argparse
from argparse import RawTextHelpFormatter
import ast
from Bio import AlignIO, SeqIO
import glob
import jenkspy
from kneebow.rotor import Rotor
import multiprocessing
import numpy as np
import os
import pandas as pd
import polars as pl
import polars.selectors as cs
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from statistics import mean


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


def all_mutations_founder(msa_clades_path, usr_metadata, path_out, ref_path):
    print("ENTERED ##########")

    meta = pl.read_csv(usr_metadata, separator="\t", quote_char=None, ignore_errors=True).fill_null(value="Undefined")
    print("METADATA OPENED!")

    all_files_msa = [msa_clades_path]

    # We are working with the sequences that have been aligned VS the reference sequence
    for ff in all_files_msa:
        clade = ff.split("/")[-1].replace(".fasta", "")

        # SEQUENCES PARSING AND CLEANING ##########

        ref_seq = SeqIO.read(ref_path, "fasta")
        ref_id = ref_seq.id

        alg = SeqIO.index(ff, "fasta")
        tot_seqs = len(alg)
        print(tot_seqs)

        # The reference sequence must be the first one of the msa
        try:
            reference = alg[ref_id]
            print(f"REFERENCE IN FASTA FILE: {reference.id} ####################")
        except KeyError:
            reference = ref_seq.seq

        headers = list(meta.get_column("Headers"))
        print(len(headers))

        # Dataframe format: {genomic position: {reference nucleotide: "", "A": n, "T": n, "C": n, "G": n}} with "n" the
        # percentage of sequences that are characterized by that mutation.
        haplos_dict = {}
        inspected_seqs = 0
        for header in headers:
            try:
                sec = alg[header]
                inspected_seqs += 1
                for i in range(0, len(reference.seq)):

                    if reference.seq[i] != sec.seq[i]:
                        pos = i+1
                        mut = sec.seq[i]

                        if pos not in haplos_dict.keys():
                            haplos_dict[pos] = {"ref": "", "A": 0, "T": 0, "C": 0, "G": 0, "-": 0}
                            haplos_dict[pos]["ref"] = reference.seq[i]
                            haplos_dict[pos][mut] += 1/tot_seqs
                        else:
                            haplos_dict[pos][mut] += 1/tot_seqs
            except KeyError:
                print(f"Header {header} not found in FASTA file")

        df_clade = pd.DataFrame.from_dict(haplos_dict).T
        df_clade.index.name = "Genomic_position"
        print(df_clade)
        print(len(df_clade))
        print(f"INSPECTED SEQS: {inspected_seqs}")
        df_clade.to_csv(f"{path_out}CL_squared_step_05/{clade}_mutation_prevalences.tsv", sep="\t")

    return


# PERCENTAGES CLUSTERING USING JENKS NATURAL BREAKS ----------
def clustering(head, tail, threshold, lower_included, clusters, pathogen_portions_file, camp, path_out):

    clade = camp.split("/")[-1].replace(".fasta", "")

    pl.Config.set_tbl_formatting("ASCII_FULL_CONDENSED")

    #  INPUTS CONTROLLER ----------
    if head:
        start_head = head-1
    else:
        start_head = 0

    if tail:
        stop_tail = tail
    else:
        stop_tail = 2000000000     # symbolic number

    if threshold:
        thr = threshold
    else:
        thr = 0.1

    if lower_included not in ["yes", "no"]:
        raise ValueError("You have entered the wrong input. User-allowed inputs: yes, no")

    # Mutation prevalence datasets
    clades_prevalences = [i for i in glob.glob(f"{path_out}CL_squared_step_05/{clade}_mutation_prevalences.tsv")]

    with open(f"{path_out}CL_squared_step_05/{clade}_cleaning_summary.txt", "w") as summary:

        # Jenks Natural Breaks breaks dictionary: {variant: [breaks]}
        per_variant_breaks = {}
        max_breaks = 0

        for cp in clades_prevalences:

            df = pl.read_csv(cp, separator="\t", ignore_errors=True).fill_null(value=0)

            summary.write("{}\nTot mutations: {}\n".format(clade, len(df)))

            indexes = df.get_column("Genomic_position")
            indexes_list = list(indexes)

            indexes = [int(i) for i in indexes_list]

            # Values to cluster: mutations percentages
            values_list = []
            discarded = 0
            for idx in sorted(indexes):    # Dataframe indexes are sequence GENOMIC POSITION)

                # If required, the user can work without sequence heads and tails
                # (recommended: less computationally demanding)
                if idx >= start_head and idx <= stop_tail:
                    temp = []
                    for col in df.columns[2:]:
                        value = df.filter(pl.col("Genomic_position") == idx).select(col).item()
                        temp.append(float(value))

                    # At now we are working with the max value
                    if max(temp)*100 >= thr:
                        # Max value of each row
                        values_list.append(max(temp))
                    else:
                        # The method is computationally demanding and we have a lot of very low values.
                        # Consequently, such nosy values are discarded
                        discarded += 1
                else:
                    discarded += 1

            summary.write("Out of range or rare mutations: {}\n".format(discarded))
            summary.write("Mutations to be clustered: {}\n".format(len(values_list)))

            if len(df)-(len(values_list) + discarded) != 0:
                raise ValueError("ERROR IN MUTATION COUNT")

            # JENKS NATURAL BREAKS IMPLEMENTATION #########################################################

            # SDAM: sum of squared deviations for array mean
            SDAM = np.sum((values_list-np.mean(values_list)) ** 2)

            # Test for different number of clusters
            GVFs = {}
            Breaks = {}     # {n_clusters: clustering breaks}
            flag = 0

            # print(len(values_list))
            # values_list = [x for x in values_list if isinstance(x, (int, float)) and not np.isnan(x)]
            # print(f"Min: {min(values_list)}, Max: {max(values_list)}, Length: {len(values_list)}")

            print(f"TOT VALUES: {len(values_list)}")
            max_clus = len(values_list)/2

            if max_clus >= 100:
                max_clus = 100

            for k in range(3, max_clus):   # len(values_list)):   # The minimum number of tested clusters is 3

                if len(values_list) <= k:
                    raise ValueError(f"Errore: Lista troppo corta. Valori: {len(values_list)}, Clusters richiesti: {k}")

                if flag == 0:

                    # Explicative example: if i have 3 ranges, the function will return 4 values:
                    #   - lower bound of the 1st class
                    #   - upper bound of the 1st class
                    #   - upper bound of the 2nd class
                    #   - upper bound of the 3rd class
                    # The algorithm evaluates autonomously the best data arrangement with k clusters
                    breaks = jenkspy.jenks_breaks(values_list, k)
                    Breaks[k] = breaks

                    # Data classification in breaks:
                    ranges = {}     # {range: [data that fall into the analysed range]}
                    for br in range(1, len(breaks)):
                        ranges["{}-{}".format(breaks[br-1], breaks[br])] = []

                    counter = 0
                    for val in values_list:
                        for rn in ranges.keys():

                            # Ranges: [A, B[
                            if lower_included == "yes":
                                if rn != list(ranges.keys())[-1]:
                                    # The lower value is included in the analysed interval while the upper value is not
                                    if val >= ast.literal_eval(rn.split("-")[0]) and val < ast.literal_eval(rn.split("-")[1]):
                                        ranges[rn].append(val)
                                        counter += 1
                                        break
                                else:
                                    # Except for the case of the last range
                                    if val >= ast.literal_eval(rn.split("-")[0]) and val <= ast.literal_eval(rn.split("-")[1]):
                                        ranges[rn].append(val)
                                        counter += 1
                                        break

                            # Ranges ]A, B]
                            elif lower_included == "no":
                                if rn != list(ranges.keys())[0]:
                                    # The upper value is included in the analysed interval while the lower value is not
                                    if val > ast.literal_eval(rn.split("-")[0]) and val <= ast.literal_eval(rn.split("-")[1]):
                                        ranges[rn].append(val)
                                        counter += 1
                                        break
                                else:
                                    # Except for the case of the first range
                                    if val >= ast.literal_eval(rn.split("-")[0]) and val <= ast.literal_eval(rn.split("-")[1]):
                                        ranges[rn].append(val)
                                        counter += 1
                                        break

                    if k == 3:
                        ref_counter = counter
                        summary.write("Classified values: {}".format(ref_counter))

                    if counter != len(values_list):
                        raise ValueError("ERROR ON CLASSIFIED VALUES NUMBER: {}".format(counter))

                    if counter != ref_counter:
                        raise ValueError("ERROR ON CLASSIFIED VALUES NUMBER: {}".format(counter))

                    # SDCM: sum of squared deviations.
                    # It is associated to the best grouping solution for k groups and, specifically, it is the
                    # lowest one.
                    SDCM = 0
                    for rn in ranges.keys():
                        for value in ranges[rn]:
                            SDCM += (value-mean(ranges[rn]))**2

                    GVFs[k] = {}
                    GVFs[k]["Value"] = (SDAM-SDCM)/SDAM

                    # Difference from previous GVF value is calculated
                    if k >= 4:
                        delta_new = GVFs[k]["Value"] - GVFs[k-1]["Value"]
                        GVFs[k]["Delta"] = delta_new

                    # When GVF value equals 1, from that moment on it will always equal 1 --> carrying capacity
                    if GVFs[k]["Value"] == 1:
                        flag = 1
                        for i in range(k+1, len(values_list)):
                            GVFs[i] = {}
                            GVFs[i]["Value"] = 1
                            GVFs[i]["Delta"] = 0
                            Breaks[i] = breaks
                else:
                    break

            GVFs[3]["Delta"] = GVFs[4]["Delta"]

            df_GVF = pd.DataFrame.from_dict(GVFs).T
            df_GVF.to_csv(f"{path_out}CL_squared_step_05/Jenks_Natural_Breaks/GVFs_{clade}.tsv", sep="\t")

            # We want to find the best number of clusters --> it is represented by the elbow of our GVF curve
            curve_data = df_GVF[['Value']]

            list_of_lists = []
            for idx in curve_data.index:
                list_of_lists.append([idx, curve_data.at[idx, 'Value']])
            array_of_arrays = np.array(list_of_lists)

            if not clusters:
                # Geometric rotation of the curve
                rotor = Rotor()
                rotor.fit_rotate(array_of_arrays)
                # Elbow point
                elbow_index = int(rotor.get_elbow_index())
                data_list = list(array_of_arrays[elbow_index])
            else:
                data_list = list(array_of_arrays[clusters])

            final_breaks = Breaks[data_list[0]]
            if len(final_breaks) > max_breaks:
                max_breaks = len(final_breaks)
            per_variant_breaks[clade] = final_breaks

            # Results plot
            fig = make_subplots(rows=1, cols=2, subplot_titles=("GVFs values", "GVF Deltas"))

            fig.add_trace(go.Scatter(x=df_GVF.index, y=df_GVF['Value'], mode="markers+lines"), row=1, col=1)
            fig.add_vline(x=data_list[0], line_width=3, line_dash="dash", line_color="green", row=1, col=1)
            fig.add_trace(go.Scatter(x=df_GVF.index, y=df_GVF['Delta'], mode="markers+lines"), row=1, col=2)

            fig.update_xaxes(title_text="Number of clusters", row=1, col=1)
            fig.update_xaxes(title_text="Number of clusters", row=1, col=2)

            fig.update_yaxes(title_text="GVF values", row=1, col=1)
            fig.update_yaxes(title_text="GVF delta values", row=1, col=2)

            fig.update_layout(title="Goodness of Variance Fit for different number of clusters and related Delta values")
            fig.write_html(f"{path_out}CL_squared_step_05/Jenks_Natural_Breaks/GVFs_plot_{clade}.html")

        for key in per_variant_breaks.keys():
            if len(per_variant_breaks[key]) < max_breaks:
                temp_list = per_variant_breaks[key]
                for i in range(len(per_variant_breaks[key]), max_breaks):
                    temp_list.append("")
                per_variant_breaks[key] = temp_list

        df_breaks = pd.DataFrame.from_dict(per_variant_breaks).T
        df_breaks.to_csv(f"{path_out}CL_squared_step_05/Jenks_Natural_Breaks/Per_Variant_Breaks_{clade}.tsv", sep="\t")

        # Clustering results plotting ----------
        dffiles = [i for i in glob.glob(f"{path_out}CL_squared_step_05/{clade}_mutation_prevalences.tsv")]

        if not os.path.isdir(f"{path_out}CL_squared_step_05/Clustered_prevalences_plots"):
            os.makedirs(f"{path_out}CL_squared_step_05/Clustered_prevalences_plots")

        df_ranges = pd.read_csv(f"{path_out}CL_squared_step_05/Jenks_Natural_Breaks/Per_Variant_Breaks_{clade}.tsv",
                                sep="\t", index_col=0).fillna("nan")

        for dffile in dffiles:

            clade = dffile.split("/")[-1].replace("_mutation_prevalences.tsv", "")

            df = pd.read_csv(dffile, sep="\t", index_col="Genomic_position")

            dict_to_plot = {}
            for idx in sorted(df.index):
                if int(idx) >= start_head and int(idx) <= stop_tail:
                    row = []
                    for col in df.columns[1:]:
                        row.append(df.at[idx, col])
                    if max(row)*100 >= thr:
                        dict_to_plot[idx] = {}
                        dict_to_plot[idx]["Max value"] = max(row)*100

            df_to_plot = pd.DataFrame.from_dict(dict_to_plot).T

            fig = px.scatter(df_to_plot, x=df_to_plot.index, y="Max value",  marginal_y="histogram")
            for column in df_ranges.columns:
                value = df_ranges.at[clade, column]
                if value != "nan":
                    fig.add_hline(y=value*100, line_width=1.5, line_color="green")
            fig.update_layout(title=f"Percentage of {clade} sequences with a specific mutation. CLUSTERED data")
            fig.write_html(f"{path_out}CL_squared_step_05/Clustered_prevalences_plots/{clade}_CLUSTERED.html")
            # ##########################################################################################################

        # Data classification on the basis of clustering results ----------
        portions = portions_dict(pathogen_portions_file)

        if not os.path.isdir(f"{path_out}CL_squared_step_05/Clustered_prevalences"):
            os.makedirs(f"{path_out}CL_squared_step_05/Clustered_prevalences")

        prevalence_tsvs = [i for i in glob.glob(f"{path_out}CL_squared_step_05/{clade}_mutation_prevalences.tsv")]

        breaks_tsv = pd.read_csv(f"{path_out}CL_squared_step_05/Jenks_Natural_Breaks/"
                                 f"Per_Variant_Breaks_{clade}.tsv", sep="\t", index_col=0).fillna("nan")

        for prevalence_tsv in prevalence_tsvs:

            df = pd.read_csv(prevalence_tsv, sep="\t", index_col=0)
            df["Orf"] = ""
            df["Group"] = ""

            clade = prevalence_tsv.split("/")[-1].replace("_mutation_prevalences.tsv", "")

            max_values = df[["A", "T", "C", "G", "-"]].max(axis=1)

            limits = []
            for col in breaks_tsv.columns:
                if breaks_tsv.at[clade, col] != "nan":
                    limits.append(breaks_tsv.at[clade, col])

            ranges = {}
            for br in range(1, len(limits)):
                ranges["{}-{}".format(limits[br-1], limits[br])] = len(limits)-br

            for idx in df.index:
                if int(idx) >= start_head and int(idx) <= stop_tail:
                    major = max_values[idx]

                    if major*100 >= thr:
                        # Data classification in ranges
                        if lower_included == "yes":
                            for rn in ranges.keys():
                                if rn != list(ranges.keys())[-1]:
                                    if major >= ast.literal_eval(rn.split("-")[0]) and major < ast.literal_eval(rn.split("-")[1]):
                                        df.at[idx, "Group"] = ranges[rn]
                                        break
                                else:
                                    if major >= ast.literal_eval(rn.split("-")[0]) and major <= ast.literal_eval(rn.split("-")[1]):
                                        df.at[idx, "Group"] = ranges[rn]
                                        break
                        elif lower_included == "no":
                            for rn in ranges.keys():
                                if rn != list(ranges.keys())[0]:
                                    if major > ast.literal_eval(rn.split("-")[0]) and major <= ast.literal_eval(rn.split("-")[1]):
                                        df.at[idx, "Group"] = ranges[rn]
                                        break
                                else:
                                    if major >= ast.literal_eval(rn.split("-")[0]) and major <= ast.literal_eval(rn.split("-")[1]):
                                        df.at[idx, "Group"] = ranges[rn]
                                        break

                        # ORF assignment
                        for portion in portions.keys():
                            try:
                                if isinstance(portions[portion], str):
                                    portions[portion] = ast.literal_eval(portions[portion])
                            except (SyntaxError, ValueError) as e:
                                print(f"Error parsing {portions[portion]}: {e}")
                            if idx in range(int(portions[portion][0]), int(portions[portion][1]+1)):
                                df.at[idx, "Orf"] = portion
                                break

            df.to_csv(f"{path_out}CL_squared_step_05/Clustered_prevalences/{clade}_haplos_prevalences_CLUSTERED.tsv",
                      sep="\t")


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="The script has the objective to recognize the most shared haplotype of a "
                                             "viral variant starting from a cleaned dataset.\nIt is optimized to exploits"
                                             " the records previously aligned VS a reference sequence which must be the"
                                             " first sequence of the package) in order to have a reference pattern of "
                                             "positions.\nATTENTION: WE SUGGEST TO USE THIS SCRIPT AFTER RUNNING THE WHOLE"
                                             " PIPELINE ON YOUR DATA WITH THE MOST POSSIBLE STRINGENT FILTERS \nYou can "
                                             "custom the clustering step as you prefer.\nThe sequences can be further cleaned"
                                             " in the following steps of the cleaning pipeline on the basis of the "
                                             "clustering results.", formatter_class=RawTextHelpFormatter)

# User FASTA files
parser.add_argument("msa_clades_path", type=str, help="FASTA file path. Insert the folder containing all your"
                                                      " files.\nThe sequences must be aligned versus a reference "
                                                      "sequence\nATTENTION: IF THE CHOSEN FOLDER CONTAINS EXTRA FASTA FILES,"
                                                      " THEY WILL BE INSPECTED BY THE SCRIPT.\nSEQUENCES SHOULD BE DIVIDED "
                                                      "INTO DIFFERENT FASTA FILES ACCORDING TO A GROUPING CRITERION: e.g. THE "
                                                      "VIRAL VARIANT THEY ARE ASSOCIATED WITH.")

# User metadata file
parser.add_argument("usr_metadata", type=str, help="Metadata path.\nThe metadata have to be correctly. If you have "
                                                   "multiple metadata files, concatenate or merge them.")

# Pathogen coding regions
parser.add_argument("pathogen_portions_file", type=str, help="Path of the .csv file containing the name of the "
                                                             "pathogen genome coding regions, each one associated "
                                                             "with the corresponding genomic positions. The file MUST "
                                                             "contain two columns: the first one named \"Name\" "
                                                             "containing the name of each coding region, the second one"
                                                             " named \"Positions\" containing the genomic positions of"
                                                             " each coding region.")

# Reference sequence path
parser.add_argument("ref_path", type=str, help="Reference sequence path.\nIf your FASTA files contain "
                                                 "the reference sequence, this parameter will be used to find the "
                                                 "reference sequence into your sequence files, otherwise it will be "
                                                 "used to compute the sequence mutations with respect to this "
                                                 "reference sequence. If you aligned your sequences in an all VS all "
                                                 "fashion, please be sure to furnish the reference sequence that was "
                                                 "aligned with all your sequences.")

# Data classification in clusters modality
parser.add_argument("lower_included", type=str, help="If this parameter is set as \"yes\", the criterion to classify "
                                                     "the data into the clusters will be: lower "
                                                     "limit is INCLUDED and upper limit is excluded with the exception "
                                                     "of the last cluster.\nIf this parameter is set as \"no\", the"
                                                     " criterion to classify data into the clusters will be: lower limit"
                                                     " is EXCLUDED and upper limit is included with the exception of "
                                                     "the first cluster.")

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

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

# Genomic position to be skipped in the head of the sequences
parser.add_argument("--head", type=int, help="Number of the genomic positions that you want to skip in the head of the "
                                             "sequences while searching for the mutations that characterize the "
                                             "analysed variant (optional):\ndefault value = 0.")

# Genomic position to be skipped in the tail of the sequences
parser.add_argument("--tail", type=int, help="Number of the genomic positions that you want to skip in the tail of the "
                                             "sequences while searching for the mutations that characterize the "
                                             "analysed variant (optional):\ndefault value = no position ignored.")

# Mutation frequency
parser.add_argument("--threshold", type=float, help="Frequence of a mutation to be considered in the clustering "
                                                    "analysis (optional).\nDefault value: 0.001 (0.1%%).")

# Number of clusters
parser.add_argument("--clusters", type=str, help="Desired number of clusters to classify your data (optional)\n."
                                                 "If you decide to NOT set this parameter, the elbow point of the "
                                                 "Goodness of Variance Fit curve will be calculated and chosen as the "
                                                 "optimal number of clusters for the classification of your data. "
                                                 "A maximum of tot_data_to_cluster/2 clusters will be tested. If this "
                                                 "value is higher than 100, a maximum of 100 clusters will be "
                                                 "evaluated.")

args = parser.parse_args()


if __name__ == "__main__":

    # INPUT CHECKING --------------------
    # FASTA files ----------
    if args.msa_clades_path.endswith(".fasta"):
        camps = [args.msa_clades_path]
    elif args.msa_clades_path.endswith("/"):
        camps = [i for i in glob.glob(f"{args.msa_clades_path}*.fasta")]
    else:
        camps = [i for i in glob.glob(f"{args.msa_clades_path}/*.fasta")]

    # Out path ----------
    if args.outfolder:
        if args.outfolder.endswith("/"):
            path_out = args.outfolder
        else:
            path_out = f"{args.outfolder}/"
    else:
        path_out = ""

    if not os.path.isdir(f"{path_out}CL_squared_step_05/"):
        os.makedirs(f"{path_out}CL_squared_step_05/")

    if not os.path.isdir(f"{path_out}CL_squared_step_05/Jenks_Natural_Breaks"):
        os.makedirs(f"{path_out}CL_squared_step_05/Jenks_Natural_Breaks")

    meta = pl.read_csv(args.usr_metadata, separator="\t", quote_char=None, ignore_errors=True).fill_null(value="Undefined")
    print("METADATA OPENED!")

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

        with open(camps[0], "r") as temp:
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
                                    + pl.col("Accession ID") + "|"
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

            if not args.index_column:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the index column")

            if not args.header_format:
                raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the header format")

            header_elements = args.header_format.split(",")

            inspected = 0

            all_idxs = list(meta.get_column(args.index_column))

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

    all_sequences = list(meta.get_column("Headers"))

    # MULTIPROCESSING
    multiprocessing.set_start_method('spawn')

    if len(camps) == 1:
        all_mutations_founder(camps[0], args.usr_metadata, path_out, args.ref_path)

    else:
        processes = []
        for camp in camps:

            p = multiprocessing.Process(target=all_mutations_founder, args=(camp,
                                                                            args.usr_metadata,
                                                                            path_out,
                                                                            args.ref_path))
            processes.append(p)

        for p in processes:
            p.start()

        for process in processes:
            process.join()

    # Multiprocessing for clustering
    if len(camps) == 1:
        clustering(args.head, args.tail, args.threshold, args.lower_included, args.clusters,
                   args.pathogen_portions_file, camps[0], path_out)

    else:
        processes = []

        for camp in camps:

            p = multiprocessing.Process(target=clustering, args=(args.head,
                                                                 args.tail,
                                                                 args.threshold,
                                                                 args.lower_included,
                                                                 args.clusters,
                                                                 args.pathogen_portions_file,
                                                                 camp,
                                                                 path_out))

            processes.append(p)

        for p in processes:
            p.start()

        for process in processes:
            process.join()
