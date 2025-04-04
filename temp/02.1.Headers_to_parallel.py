# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - gisaid type: if you downloaded your sequences from the GISAID database you must specify if you have downloaded
#                them from the "Downloads" section (A) or from the "Search" section (B).
#		 ATTENTION: if you downloaded your sequences from the "Downloads" section please download also your
#		            metadata from the "Downloads" section.
# - metadata:
#	mandatory columns: Accession ID, Collection Date, Submission date, Virus name (GISAID)
#                          Accession (NCBI)
#			   Columns containing the data contained in the sequence header (Personal metadata).
#                          See "header format"
#
# - number of threads: number of available threads to analyze your data. The total number of sequences will be
#                      split into a number of FASTA files equal to the number of available threads in order to implement
#                      a multiprocessing pipeline.
#
# - header format: If your metadata have not GISAID or Ncbi format, you can use them anyway but you must indicate the
#                  format of the headers of your sequences. If your header has the following format:
#
#                  	"Strain/Accession_id/Collection_date"
#
#                  the --header_format string MUST be like this:
#
#                  	"\"STRAIN COLUMN NAME\", \"/\" (SEPARATOR), \"ACCESSION ID COLUMN NAME\", \"/\", \"COLLECTION DATE COLUMN NAME\"
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
#                  Please pay attention at the white spaces. If your headers contain white spaces, please remember that
#		   all the information contained in the header after the first white space will be considered as part
#		   of the sequence description when parsed. It is preferable to remove or substitute all teh white spaces.
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# This script exploits sequences metadata to generate .txt files containing the reconstructed headers of the FASTA
# file. The number of output files corresponds to the number of threads that you want to use and must be jointly used
# with the "parallel_faidx.sh" script to generate the corresponding FASTA files starting from a single big file.
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - .txt files containing the headers to obtain a faster multiprocessing pipeline.
# - Summary (.txt) file.
# ----------------------------------------------------------------------------------------------------------------------
#
# ######################################################################################################################


import argparse
from argparse import RawTextHelpFormatter
import math
import os
import pandas as pd
import polars as pl


# FUNCTIONS ############################################################################################################
def header_formatter(usr_metadata, threads, gisaid_type, index_column, header_format, outfolder):

	#  INPUTS CONTROLLER ----------
	if outfolder:
		if outfolder.endswith("/"):
			path_out = outfolder
		else:
			path_out = f"{outfolder}/"
	else:
		path_out = ""

	if not os.path.isdir(f"{path_out}CL_squared_parallel"):
		os.makedirs(f"{path_out}CL_squared_parallel")

	all_seqs = []
	with open(f"{path_out}CL_squared_parallel/Summary.txt", "w") as summary:

		cleaned_meta = pl.read_csv(usr_metadata, separator="\t", ignore_errors=True).fill_null(value="Undefined")
		summary.write(f"Metadata length: {len(cleaned_meta)}\n")

		meta_columns = cleaned_meta.columns

		number_of_sec_per_file = math.ceil(len(cleaned_meta) / threads)
		summary.write(f"Number of sequences per file (threads based): {number_of_sec_per_file}\n")

		type = ""

		# Public databanks metadata
		if not header_format:
			if "Accession ID" not in meta_columns and "Accession" not in meta_columns:
				raise ValueError("Wrong column names or missing columns. Minimum mandatory "
						 "columns for GISAID or NCBI metadata: see info at --help")

			elif "Accession ID" in meta_columns:
				type = "Gisaid"
				print(type)

			elif "Accession" in meta_columns:
				type = "Ncbi"
				print(type)

		# Private metadata
		else:
			type = "personal"
			print(type)

		if type == "Gisaid":
			if not gisaid_type:
				raise ValueError("Specify the type of the data you downloaded.\n"
						 "Admitted ones: \"A\" for comprehensive FASTA files retrieved from DOWNLOAD SECTION\n"
						 "\"B\" for FASTA files retrieved from SEARCH SECTION\nATTENTION: if you"
						 " download your sequences from this section, please download ALSO YOUR "
						 "METADATA FROM THIS SECTION")

			elif gisaid_type not in ["A", "B"]:
				raise ValueError("Wrong type. "
						 "Admitted ones: \"A\" for comprehensive FASTA files retrieved from DOWNLOAD SECTION\n"
						 "\"B\" for FASTA files retrieved from SEARCH SECTION. ATTENTION: if you"
						 " download your sequences from this section, please download ALSO YOUR "
						 "METADATA FROM THIS SECTION")

			if gisaid_type == "A" and "Headers" not in meta.columns:

				cleaned_meta = cleaned_meta.with_columns((pl.col("Virus name") + "|"
                                                                        + pl.col("Collection date")+ "|"
                                                                        + pl.col("Submission date")).alias('Headers'))

				cleaned_meta = cleaned_meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
				cleaned_meta = cleaned_meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
				cleaned_meta = cleaned_meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))

			elif gisaid_type == "B" and "Headers" not in meta.columns:


				cleaned_meta = cleaned_meta.with_columns((pl.col("Virus name") + "|"
                                                                        + pl.col("Accession ID")+ "|"
                                                                        + pl.col("Collection date")).alias('Headers'))

				cleaned_meta = cleaned_meta.with_columns(pl.col("Headers").str.replace_all("|Undefined", ""))
				cleaned_meta = cleaned_meta.with_columns(pl.col("Headers").str.replace_all("Undefined", ""))
				cleaned_meta = cleaned_meta.with_columns(pl.col("Headers").str.replace_all(" ", "_"))

			inspected = len(cleaned_meta)

			all_seqs = list(cleaned_meta.get_column("Headers"))

			if len(all_seqs) != len(cleaned_meta):
				summary.write("ATTENTION: DUPLICATE HEADERS")
				print("ATTENTION: DUPLICATE HEADERS")

		elif type == "Ncbi":

			cleaned_meta = cleaned_meta.with_column(pl.col("Accession").str.replace_all(" ", "_"))

			all_seqs = list(cleaned_meta.get_column("Accession"))

			if len(all_seqs) != len(cleaned_meta):
				summary.write("ATTENTION: DUPLICATE HEADERS")
				print("ATTENTION: DUPLICATE HEADERS")

			inspected = len(cleaned_meta)

		else:
			if not index_column:
				raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the index column")

			if not header_format:
				raise ValueError("If your metadata have not GISAID or Ncbi format, you must indicate the header format")

			header_elements = header_format.split(",")

			inspected = 0

			all_idxs = list(cleaned_meta.get_column(index_column))

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

		summary.write(f"Inspected sequences: {inspected}")

		# Headers files
		print(len(all_seqs))
		j = 0
		for i in range(0, len(all_seqs), number_of_sec_per_file):
			group = all_seqs[i:i + number_of_sec_per_file]
			j += 1

			if not os.path.isdir(f"{path_out}CL_squared_parallel/Headers_Groups"):
				os.makedirs(f"{path_out}CL_squared_parallel/Headers_Groups")

			with open(f"{path_out}CL_squared_parallel/Headers_Groups/Headers_group_{j}.txt", "w") as doc:
				for el in group:
					doc.write(f"{el}\n")

		cleaned_meta.write_csv(f"{path_out}CL_squared_parallel/Cleaned_meta_with_headers.tsv", separator="\t")


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="This script exploits sequences metadata to generate .txt files containing"
                                             " the reconstructed headers of the FASTA file. The number of output files "
                                             "corresponds to the number of threads that you want to use and must be "
                                             "jointly used with the \"parallel_faidx.sh\" script to generate the "
                                             "corresponding FASTA files starting from a single big file.",
                                             formatter_class=RawTextHelpFormatter)

# User metadata file
parser.add_argument("usr_metadata", type=str, help="Metadata file path.\nThe metadata have to correctly be "
                                                   "associated with the sequences that you want to inspect.\nIf you "
                                                   "have multiple metadata files, concatenate or merge them.\nMANDATORY FORMAT: .tsv\n"
                                                   "Mandatory columns:\nAccession ID, Collection Date, Submission date, Virus name "
                                                   "(GISAID)\nAccession (NCBI)\nColumns containing the data contained in the sequence "
                                                   "header (Personal metadata). See --header format.\nATTENTION: When downloading "
                                                   "METADATA from NCBI, PLEASE DOWNLOAD THE METADATA CONTAINING THE GENOMIC SEQUENCE"
                                                   " VERSION. MOREOVER, WHEN YOU DOWNLOAD THE SEQUENCES FILE PLEASE UNPIN ALL THE "
                                                   "COLUMNS DIFFERENT FROM \"ACCESSION\" (the header must be without description)")

# NUmber of threads
parser.add_argument("threads", type=int, help="Number of threads that you want to use in the next steps of the pipeline.")

# Gisaid downloaded data type
parser.add_argument("--gisaid_type", type=str, help="Type of the FASTA FILE you downloaded from GISAID database. "
                                                    "MANDATORY IF YOUR SEQUENCES COME FROM THIS REPOSITORY.\n"
                                                    "Admitted ones: \"A\" for comprehensive FASTA files retrieved from "
                                                    "DOWNLOAD SECTION\nATTENTION: if you download your "
                                                    "sequences from this section, please download ALSO YOUR "
                                                    "METADATA FROM THIS SECTION.\n\"B\" for FASTA files retrieved from SEARCH SECTION.")

# Index column for metadata different from GISAID, Ncbi
parser.add_argument("--index_column", type=str, help="If your metadata have not GISAID or Ncbi format, you can use them"
                                                     " anyway but you must indicate the index column of the file.")

# Header format for metadata different from GISAID, Ncbi
parser.add_argument("--header_format", type=str, help="If your metadata have not GISAID or Ncbi format, you can use "
                                                      "them anyway but you must indicate the format of the headers of "
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
                                                      " spaces. Please pay attention at the white spaces. If your headers "
                                                      "contain white spaces, please remember that all the information "
                                                      "contained in the header after the first white space will be "
                                                      "considered as part of the sequence description when parsed. It is "
                                                      "preferable to remove or substitute all teh white spaces.\nATTENTION: "
                                                      "This option could be slower.")

# Outputs path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

args = parser.parse_args()

header_formatter(args.usr_metadata, args.threads, args.gisaid_type, args.index_column, args.header_format, args.outfolder)
