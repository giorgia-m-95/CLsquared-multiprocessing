# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - Metadata file.
#   MANDATORY MINIMUM FORMAT:
#    - columns: Accession ID (sequence ID), Collection date, Location, Host (GISAID metadata)
#               Accession (sequence ID), Collection_Date, Country, Host (NCBI metadata - please download the accession
#               ID WITH VERSION)
#
#               If your metadata file doesn't match the GISAID and NCBI formats, you MUST provide
#               a list of columns as an extra input. The list should follow this pattern:
#
#               ["sequence id column name", "host column name", "place column name", "collection date column name"].
#
#               ATTENTION: if there is a missing a column, do maintain this order and use "" in its place.
#               For instance, if you lack the "host" column, your list should look like this:
#               ["sequence id column name", "", "place column name", "collection date column name"].
#
#               It is MANDATORY TO HAVE A SEQUENCE ID COLUMN.
#
# - Completeness
#   The user can decide wether to work only with complete collection dates or to retain the incomplete collection dates
#   that will be properly filled.
#
# - Out folder to save the summary file and the cleaned metadata (optional). Default out folder: current folder.
#
# - Place (optional).
#   The user can decide to work only with the sequences collected in a specific place.
#
# - Start date/end date (optional).
#   The user can decide to filter the sequences on the basis of the collection date.
#   MANDATORY DATE FORMATS: "YYYY-MM-DD", "YYYY-MM", "YYYY".
#
# - Host (optional)
#   The user can decide to filter the sequences on the basis of the host.
#
# IN ORDER TO BETTER UNDERSTAND HOW TO CHOOSE THE CLEANING PARAMETERS RUN THE
# 01.Metadata_cleaner_PRE_CLEANING_DIAGNOSTICS_WITH_OLD_VERSION.py SCRIPT AND READ THE PRE CLEANING DIAGNOSTICS.
#
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# Substitution of non-ASCII encoded characters with ASCII encoded ones in the location column.
#
# IF REQUIRED by the user:
#  - Selection of "human" or any other desired host sequences.
#  - Selection of the sequences associated with a complete collection date.
#  - Selection of the sequences sampled in a required place.
#  - Selection of the sequences sampled in a required time span.
#  - Writing of the .txt files containing the Accession IDs of the selected sequences.
#    If you are working with the GISAID database, this file format is suitable for the download of the sequences from
#    the GISAID website by simply uploading it (max 10.000 sequences, as required by the policies of GISAID).
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - Cleaned metadata (.tsv) file
# - Accession IDs files (if required by the user)
# - Summary (.txt) file
# ----------------------------------------------------------------------------------------------------------------------
#
# ######################################################################################################################


import argparse
from argparse import RawTextHelpFormatter
import ast
from datetime import datetime
import os
import polars as pl


# FUNCTIONS ############################################################################################################

def meta_cleaner(usr_metadata, outfolder, host, place, start_date, end_date, completeness, meta_columns_list,
                 duplicates):

    pl.Config.set_tbl_formatting("ASCII_FULL_CONDENSED")

    #  INPUTS CONTROLLER ----------
    if outfolder:
        if outfolder.endswith("/"):
            path_out = outfolder
        else:
            path_out = f"{outfolder}/"
    else:
        path_out = ""

    if completeness not in ["y", "n"]:
        raise ValueError("Wrong input. User-allowed inputs: y, n")

    if duplicates and duplicates not in ["y", "n"]:
            raise ValueError("Wrong input. User-allowed inputs: y, n")

    if start_date:
        if not len(start_date) in [4, 7, 10]:
            raise ValueError("Wrong date format. Mandatory formats: \"YYYY-MM-DD\", \"YYYY-MM\" or \"YYYY\"")
        else:
            if len(start_date) == 4:
                start_date += "-01-01"
            if len(start_date) == 7:
                start_date += "-15"
    if end_date:
        if not len(end_date) in [4, 7, 10]:
            raise ValueError("Wrong date format. Mandatory format: \"YYYY-MM-DD\", \"YYYY-MM\" or \"YYYY\"")
        else:
            if len(end_date) == 4:
                end_date += "-01-01"
            if len(end_date) == 7:
                end_date += "-15"

    # Out folder creation ----------
    if not os.path.isdir(f"{path_out}CL_squared_step_01"):
        os.makedirs(f"{path_out}CL_squared_step_01")

    # Metadata filtering ---------
    with open(f"{path_out}CL_squared_step_01/Metadata_filtering_summary.txt", "w") as doc:

        meta = pl.read_csv(usr_metadata, separator="\t", quote_char=None, ignore_errors=True).fill_null(value="Undefined")
        print("METADATA OPENED: STARTING DIAGNOSTICS...")

        meta_columns = list(meta.columns)

        doc.write(f"STARTING METADATA LENGTH: {len(meta)}\n\n")

        # Input columns checking ----------
        # Public databanks metadata
        if not meta_columns_list:
            if "Accession ID" not in meta_columns and "Accession" not in meta_columns:
                raise ValueError(
                    "Wrong column names or missing columns. Minimum mandatory columns for GISAID or NCBI metadata:"
                    " Sequence accession ID, Collection date, Location")

            elif "Accession ID" in meta_columns:
                type = "Gisaid"

            elif "Accession" in meta_columns:
                type = "Ncbi"

        # Private metadata
        else:
            type = "personal"

            meta_columns = ast.literal_eval(meta_columns_list)

            if len(meta_columns) != 4:
                raise ValueError(
                    "Missing columns! Mandatory list format: "
                    "['index column name', 'host column name', 'place column name', "
                    "'collection date column name']. If you have a missing column, do maintain the indicated order"
                    " and use '' in its place.")

        # Checking for duplicate entries if Virus name in metadata columns: VALID FOR GISAID METADATA ----------
        if type == "Gisaid":

            if "Virus name" in meta_columns:

                vir_names = meta.get_column("Virus name")
                vir_names = list(set(vir_names))

                if len(meta) != len(vir_names):

                    doc.write(f"Copies of existing strains: {len(meta)-len(vir_names)}\n")

                    if duplicates == "y":

                        meta = meta.unique(subset=["Virus name"], keep="none")
                        doc.write("Not unique records dropped\n")

                    else:
                        doc.write("Not unique records retained\n")

                else:
                    doc.write("No duplicate strains\n")

                doc.write(f"Unique strains checked. Survived unique metadata entries: {len(meta)}\n\n")

        # Removal of sequences with incomplete collection dates ----------
        doc.write("\nCollection date completeness check ----------\n")

        # If the user wants to retain sequences with missing collection date:
        if type == "Gisaid":
            cd = "Collection date"
        elif type == "Ncbi":
            if "Collection_Date" in meta_columns:
                cd = "Collection_Date"
            else:
                cd = "Isolate Collection date"
        else:
            cd = meta_columns[3]

        if cd != "":

            meta = meta.with_columns(pl.col(cd).str.replace_all("Undefined", ""))
            meta = meta.with_columns(pl.col(cd).str.replace_all("/", "-"))
            meta = meta.with_columns(pl.col(cd).str.replace_all("-XX", ""))

            meta_not_complete_date = meta.filter(pl.col(cd).str.len_chars() != len("YYYY-MM-DD"))

            if len(meta_not_complete_date) != 0:

                doc.write(f"Records associated with an incomplete collection date: {len(meta_not_complete_date)}\n")

                if completeness == "y":

                    meta = meta.filter(pl.col(cd).str.len_chars() == len("YYYY-MM-DD"))

                    doc.write("Records associated with an incomplete collection date: dropped\n\n")

                else:

                    doc.write("Records associated with an incomplete collection date: retained\n\n")

            else:

                doc.write("No incomplete collection dates\n\n")

            doc.write(f"Surviving date-filtered metadata (based on completeness): {len(meta)}\n\n")

        # Host selection ----------
        if host:
            if type == "Gisaid":
                host_col = "Host"
            elif type == "Ncbi":
                if "Host" in meta_columns:
                    host_col = "Host"
                else:
                    host_col = "Host Name"
            else:
                host_col = meta_columns[1]

            host_name = host.title()

            hosts = meta.get_column(host_col)
            all_hosts = hosts.str.to_titlecase()

            hosts_list = list(set(all_hosts))

            if host_name in hosts_list:
                meta = meta.filter(pl.col(host_col).str.to_titlecase() == host_name)
            else:
                raise ValueError(f"The host name you entered is not available in your metadata. Available hosts: "
                                 f"{list(set(hosts_list))}")

            doc.write(f"\nSurvived metadata according to host selection (host = {host_name}): {len(meta)}\n\n")

        # Editing of country names ----------
        if type == "Gisaid":
            loc = "Location"
        elif type == "Ncbi":
            if "Country" in meta_columns:
                loc = "Country"
            else:
                loc = "Submitter Country"
        else:
            loc = meta_columns[2]

        if loc != "":

            # GISAID Location format: Continent/Country/Place specifications
            if type == "Gisaid":

                locations = list(meta.get_column("Location").str.to_titlecase())

                new_locations = []
                for i in locations:
                    if len(i.split(" / ")) >= 2:
                        new_locations.append(i.split(" / ")[1])
                    else:
                        new_locations.append("Undefined")

                idx = meta_columns.index("Location")

                new_col = pl.Series("Location", new_locations)
                meta = meta.replace_column(idx, new_col)

                meta = meta.fill_null(value="Undefined")

            replacement_dict = {" ": "_", "é": "e", "è": "e", "ô": "o"}  # no ASCII encoded characters

            for key in replacement_dict:
                meta = meta.with_columns(pl.col(loc).str.replace_all(key, replacement_dict[key]))

            # Country selection ----------
            if place:

                place_name = place.title()

                places_list = list(set(meta.get_column(loc).str.to_titlecase()))

                if place_name in places_list:
                    meta = meta.filter(pl.col(loc).str.to_titlecase() == place_name)
                else:
                    similar_places = [el for el in places_list if place_name in el]

                    if len(similar_places) != 0:
                        raise ValueError(f"The place you entered is not available in your metadata. "
                                         f"COMPLETE PLACES LIST: {places_list}.\nSIMILAR PLACES: {similar_places}")
                    else:
                        raise ValueError(f"The place you entered is not available in your metadata. "
                                         f"Complete places list: {places_list}")

                doc.write(
                    f"\nSurvived metadata according to place selection (place = {place_name}): {len(meta)}\n\n")

        # Date filtering ----------
        if cd != "":

            # Correct date format

            meta = meta.filter((pl.col(cd).str.len_chars() == 4) |
                               (pl.col(cd).str.len_chars() == 7) |
                               (pl.col(cd).str.len_chars() == 10))

            doc.write(f"Sequences with wrong format collection date dropped. Filtered metadata length: {len(meta)}\n\n")

            # Dates completeness

            if completeness == "n":

                meta = meta.with_columns(pl.col(cd).map_elements(lambda x: x + "-01-01" if len(x) == 4 else x, return_dtype=str).alias("Filled_dates"))
                meta = meta.with_columns(pl.col("Filled_dates").map_elements(lambda x: x + "-15" if len(x) == 7 else x, return_dtype=str))

            elif completeness == "y":

                meta = meta.filter(pl.col(cd).str.len_chars() == len("YYYY-MM-DD"))

            # Time ranges

            if start_date and end_date:

                start = datetime.strptime(start_date, '%Y-%m-%d')
                end = datetime.strptime(end_date, '%Y-%m-%d')

                if completeness == "y":
                    meta = meta.filter(pl.col(cd) != "Undefined")
                    meta = meta.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                    meta = meta.filter((pl.col(cd) >= start) &
                                       (pl.col(cd) <= end))
                else:
                    meta = meta.filter(pl.col("Filled_dates") != "Undefined")
                    meta = meta.with_columns(pl.col("Filled_dates").str.to_date("%Y-%m-%d"))
                    meta = meta.filter((pl.col("Filled_dates") >= start) &
                                       (pl.col("Filled_dates") <= end))

                doc.write(f"\nLength of the metadata filtered on the basis of the collection date: {len(meta)}\n\n")

            elif start_date:

                start = datetime.strptime(start_date, '%Y-%m-%d')

                if completeness == "y":
                    meta = meta.filter(pl.col(cd) != "Undefined")
                    meta = meta.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                    meta = meta.filter(pl.col(cd) >= start)
                else:
                    meta = meta.filter(pl.col("Filled_dates") != "Undefined")
                    meta = meta.with_columns(pl.col("Filled_dates").str.to_date("%Y-%m-%d"))
                    meta = meta.filter(pl.col("Filled_dates") >= start)

                doc.write(f"\nLength of the metadata filtered on the basis of the collection date: {len(meta)}\n\n")

            elif end_date:

                end = datetime.strptime(end_date, '%Y-%m-%d')

                if completeness == "y":
                    meta = meta.filter(pl.col(cd) != "Undefined")
                    meta = meta.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                    meta = meta.filter(pl.col(cd) <= end)
                else:
                    meta = meta.filter(pl.col("Filled_dates") != "Undefined")
                    meta = meta.with_columns(pl.col("Filled_dates").str.to_date("%Y-%m-%d"))
                    meta = meta.filter(pl.col("Filled_dates") <= end)

                doc.write(f"\nLength of the metadata filtered on the basis of the collection date: {len(meta)}\n\n")

        # Saving meta TSV ----------
        doc.write(f"\nFINAL METADATA LENGTH {len(meta)}")

        meta.write_csv(f"{path_out}CL_squared_step_01/Cleaned_metadata.tsv", separator="\t")

        return meta, path_out, type, meta_columns_list


def epi_files_maker(meta, path_out, type, meta_columns_list):

    if type == "Gisaid":
        idx = "Accession ID"
    elif type == "Ncbi":
        idx = "Accession"
    else:
        idx = ast.literal_eval(meta_columns_list)[0]

    to_download_epis = list(meta.get_column(idx))

    with open(f"{path_out}CL_squared_step_01/Accession_IDs_files_summary.txt", "w") as doc:
        j = 0
        counter_added = 0
        for epi in to_download_epis:
            counter_added += 1

            # EPIs files for sequence download can't contain more than 10.000 EPIs
            if counter_added % 10000 != 0:
                with open(f"{path_out}CL_squared_step_01/EPIs_list_{j}.txt", "a") as epi_doc:
                    epi_doc.write(f"{epi}\n")
            else:
                with open(f"{path_out}CL_squared_step_01/EPIs_list_{j}.txt", "a") as epi_doc:
                    epi_doc.write(f"{epi}\n")

                j += 1

        doc.write(f"Metadata length: {len(meta)}\n\n")
        doc.write(f"Tot written Accession IDs: {counter_added}\n\n")
        doc.write(f"Tot written files: {j+1}\n\n")


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="Substitution of non-ASCII encoded characters with ASCII encoded ones in "
                                             "the location column.\n"
                                             "IF REQUIRED by the user:\n"
                                             "- Selection of human or any other desired host sequences.\n"
                                             "- Selection of the sequences associated with complete collection dates.\n"
                                             "- Selection of the sequences sampled in a required place.\n"
                                             "- Selection of the sequences sampled in a required time span.\n"
                                             "- Writing of the .txt files containing the accession IDs of the selected "
                                             "sequences. "
                                             "If you are working with the GISAID website, this file format is suitable"
                                             " for the download of the sequences by simply uploading it (max 10000 "
                                             "sequences, as required by the policies of GISAID).\n",
                                 formatter_class=RawTextHelpFormatter)

# User metadata file
parser.add_argument("usr_metadata", type=str, help="Metadata file path.\nSequence ID, Collection date, Location "
                                                   "and Host columns are mandatory.\nIf you work with NCBI metadata"
                                                   " please download the accession ID WITH VERSION.\nIf you have "
                                                   "multiple metadata files, concatenate or merge them. MANDATORY "
                                                   "FILE FORMAT: .tsv")

# Date completeness
parser.add_argument("completeness", type=str, help="If you want the collection dates to be complete (YYYY-MM-DD "
                                                     "format) set this parameter as \"y\", otherwise set is as "
                                                     "\"n\".")

# Output folder path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

# Virus host
parser.add_argument("--host", type=str, help="Selection of the sequences associated with the required host (optional)."
                                             "\nIf not specified, all the available hosts will be retained.")

# Place
parser.add_argument("--place", type=str, help="Selection of the sequences sampled in a required place (optional).\nIf "
                                              "not specified, all the available places will be retained.")


# Start date
parser.add_argument("--start_date", type=str, help="Selection of the sequences sampled after the required collection "
                                                   "date (optional).\nIf you want to crate a specific time span "
                                                   "add the --end_date optional parameter. Mandatory date format: "
                                                   "\"YYYY-MM-DD\".")

# End date
parser.add_argument("--end_date", type=str, help="Selection of the sequences sampled before the required collection "
                                                 "date (optional).\nIf you want to crate a specific time span "
                                                 "add the --start_date optional parameter. Mandatory date format: "
                                                 "\"YYYY-MM-DD\".")

# Metadata columns
parser.add_argument("--meta_columns_list", type=str, help="If your metadata file doesn't match the GISAID and NCBI "
                                                          "formats, you MUST provide a list of columns as an extra input"
                                                          ". The list should follow this pattern:\n[\"sequence id "
                                                          "column name\", \"host column name\", \"place column name\","
                                                          " \"collection date column name\"]. Please use double quotes "
                                                          "before and after the square brakets and signle quotes before and"
                                                          " after every list element to correctly insert the input variable."
                                                          "\n ATTENTION: if you have "
                                                          "a missing column, maintain the order and use \"\" in its "
                                                          "place. For instance, if you lack the \"host\" column, your "
                                                          "list should look like this: [\"sequence id column name\", "
                                                          "\"\", \"place column name\", \"collection date column name"
                                                          "\"]. It is MANDATORY TO HAVE A SEQUENCE ID COLUMN.")

# Viral strain duplicates
parser.add_argument("--duplicates", type=str, help="If you work with GISAID database and you want the sequences "
                                                   "associated with duplicated viral strains to be removed, set this"
                                                   " parameter as \"y\", otherwise set is as \"n\" (optional).")

args = parser.parse_args()

cleaned_meta, path_out, type, meta_columns_list = meta_cleaner(args.usr_metadata, args.outfolder, args.host, args.place,
                                                               args.start_date, args.end_date, args.completeness,
                                                               args.meta_columns_list, args.duplicates)

# If the user gives permission, txt files containing the accession IDs associated with the survived sequences are
# written (10.000 records per file).
request = input("\n\n ################################################################################################"
                " Do you want to write the .txt files containing the survived sequences accession IDs?\nThe files will"
                " have a format suitable to download the sequences associated with the IDs from GISAID website by "
                "simply uploading them (10.000 accession IDs at most per file).\n(y/n): ")

if request == "y":
    epi_files_maker(cleaned_meta, path_out, type, meta_columns_list)
elif request == "n":
    print("01.Fasta_writer.py ended without errors!")
else:
    raise ValueError("Wrong input. User-allowed inputs: y, n")
