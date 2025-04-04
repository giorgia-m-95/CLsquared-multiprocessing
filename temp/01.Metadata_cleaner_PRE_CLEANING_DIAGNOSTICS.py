# INFO #################################################################################################################
#
# INPUTS ---------------------------------------------------------------------------------------------------------------
# - Metadata file to be cleaned. MANDATORY FORMAT: tsv FILE
#
# - Out folder to save the summary file with the diagnostics file (optional).
#   Default out folder: current folder.
# ----------------------------------------------------------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------------------------------------------------------
# The script inspects the metadata file to filter in order to provide useful information to consciously clean the
# dataset.
# Specifically, information with respect to the parameters that you will choose in the cleaning step is contained in the
# output file:
#   - Host
#   - Place
#   - Sequences collected after a "start date", before an "end date" or between a "start" and an "end date"
#
# If your metadata file doesn't match the GISAID or NCBI formats, you can provide a list of columns as an extra input.
# The list must follow this pattern:
#       ["sequence id column name", "host column name", "place column name", "collection date column name"].
#
# ATTENTION: if there is a missing column, do maintain this order and use "" in its place. For example, if you lack
# the "host" column, your list should look like this:
#       ["sequence id column name", "", "place column name", "collection date column name"].
#
# ----------------------------------------------------------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------------------------------------------------------
# - Diagnostics file
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
def meta_cleaner_diagnostic(usr_metadata, outfolder, host, place, start_date, end_date, meta_columns_list):

    #  INPUTS CONTROLLER ----------
    if outfolder:
        if outfolder.endswith("/"):
            path_out = outfolder
        else:
            path_out = f"{outfolder}/"
    else:
        path_out = ""

    if start_date:
        if not len(start_date) in [4, 7, 10]:
            raise ValueError("Wrong date format. Mandatory formats: \"YYYY-MM-DD\", \"YYYY-MM\" or \"YYYY\"")
        else:
            if len(start_date) == 4:
                start_date += "-01-01"
            elif len(start_date) == 7:
                start_date += "-15"

    if end_date:
        if not len(end_date) in [4, 7, 10]:
            raise ValueError("Wrong date format. Mandatory format: \"YYYY-MM-DD\", \"YYYY-MM\" or \"YYYY\"")
        else:
            if len(end_date) == 4:
                end_date += "-01-01"
            elif len(end_date) == 7:
                end_date += "-15"

    # Out folder creation ----------
    if not os.path.isdir(f"{path_out}CL_squared_step_01_DIAGNOSTICS"):
        os.makedirs(f"{path_out}CL_squared_step_01_DIAGNOSTICS")

    # Metadata pre filtering diagnostics ---------
    with open(f"{path_out}CL_squared_step_01_DIAGNOSTICS/Metadata_pre_filtering_diagnostics.txt", "w") as doc:

        meta = pl.read_csv(usr_metadata, separator="\t", ignore_errors=True).fill_null(value="Undefined")
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
                    "['ids column name', 'host column name', 'place column name', "
                    "'collection date column name']. If you have a missing column, do maintain the indicated order"
                    " and use '' in its place.")

        # Detection of sequences with some missing information in metadata (Collection date, Location, Host) ----------
        if type == "Gisaid":

            # Checking for duplicate entries if Virus name in metadata columns.
            # The number of duplicate entries is reported.
            if "Virus name" in meta_columns:

                vir_names = meta.get_column("Virus name")
                vir_names = list(set(vir_names))

                if len(meta) != len(vir_names):
                    doc.write(f"Copies of existing strains: {len(meta) - len(list(set(vir_names)))}\n\n")

                else:
                    doc.write("No duplicate strains\n\n")

            meta_temp = meta.select(pl.col(["Collection date", "Location", "Host"]))

            meta_temp = meta_temp.filter((pl.col('Collection date') == 'Undefined') |
                                         (pl.col('Location') == 'Undefined') |
                                         (pl.col('Host') == 'Undefined'))

        elif type == "Ncbi":

            meta_temp = meta.select(pl.col(["Collection_Date", "Country", "Host"]))

            meta_temp = meta_temp.filter((pl.col('Collection_Date') == 'Undefined') |
                                         (pl.col('Country') == 'Undefined') |
                                         (pl.col('Host') == 'Undefined'))

        else:

            temp_meta_cols = [i for i in meta_columns if i != ""]

            meta_temp = meta.select(pl.col(temp_meta_cols))

            indexes_to_cons = []
            for col in temp_meta_cols[1:]:
                mini_meta = meta_temp.filter(pl.col(col) == 'Undefined')
                indexes_to_cons.extend(list(mini_meta.get_column(meta_columns[0])))
                indexes_to_cons = list(set(indexes_to_cons))

            meta_temp = meta_temp.filter(pl.col(meta_columns[0]).is_in(indexes_to_cons))

        missing_meta_seqs = len(meta_temp)

        if missing_meta_seqs != 0:
            doc.write(f"Sequences with MISSING Collection date, Location and/or Host: {missing_meta_seqs}\n\n")
        else:
            doc.write(f"No sequences with MISSING Collection date, Location and/or Host\n\n")

        # Incomplete collection dates detection ----------
        if type == "Gisaid":
            cd = "Collection date"
        elif type == "Ncbi":
            cd = "Collection_Date"
        else:
            cd = meta_columns[3]

        if cd != "":

            meta_temp = meta.select(meta_columns)

            dates_temp = pl.Series(cd, meta_temp.get_column(cd))

            dates_temp = dates_temp.str.replace("Undefined", "")
            dates_temp = dates_temp.str.replace_all("/", "-")
            dates_temp = dates_temp.str.replace_all("-XX", "")

            dates_temp = list(dates_temp)
            dates_temp = [i for i in dates_temp if len(i) != len("YYYY-MM-DD")]

            if len(dates_temp) != 0:
                doc.write(f"Records associated with an incomplete collection date: {len(dates_temp)}\n\n")
            else:
                doc.write("No incomplete collection dates\n\n")

        # Host selection ----------
        if type == "Gisaid":
            host_col = "Host"
        elif type == "Ncbi":
            host_col = "Host"
        else:
            host_col = meta_columns[1]

        if host_col != "":

            if host:

                host_name = host.title()

                all_hosts = list(meta.get_column(host_col).str.to_titlecase())
                hosts_list = list(set(all_hosts))

                if host_name in hosts_list:
                        selected_host = [i for i in all_hosts if i == host_name]
                else:
                    raise ValueError(f"The host name you entered is not available in your metadata. "
                                     f"Available hosts: {hosts_list}")

                doc.write(f"Sequences associated with selected host (host = {host_name}): {len(selected_host)}\n\n")

        # Country selection ----------
        if type == "Gisaid":
            loc = "Location"
        elif type == "Ncbi":
            loc = "Country"
        else:
            loc = meta_columns[2]

        if loc != "":

            if place:

                place_name = place.title()

                locations = list(meta.get_column(loc).str.to_titlecase())

                # GISAID Location format: Continent/Country/Place specifications
                if type == "Gisaid":
                    locations = [i.split(" / ")[1] for i in locations]

                places_list = list(set(locations))

                if place_name in places_list:
                    selected_location = [i for i in locations if i == place_name]
                else:
                    similar_places = [el for el in places_list if place_name in el]

                    if len(similar_places) != 0:
                        raise ValueError(f"The place you entered is not available in your metadata. "
                                         f"COMPLETE PLACES LIST: {places_list}.\nSIMILAR PLACES: {similar_places}")
                    else:
                        raise ValueError(f"The place you entered is not available in your metadata. "
                                         f"Complete places list: {places_list}")
                doc.write(
                    f"Sequences associated with selected place (place = {place_name}): {len(selected_location)}\n\n")

        # Date filtering ----------
        if cd != "":

            # Dates filling
            if start_date and end_date:

                # Dates formatting
                start = datetime.strptime(start_date, '%Y-%m-%d')
                end = datetime.strptime(end_date, '%Y-%m-%d')

                temp_meta = meta.filter(pl.col(cd) != "Undefined")

                temp_meta = temp_meta.with_columns(pl.col(cd).str.replace_all("/", "-"))
                temp_meta = temp_meta.with_columns(pl.col(cd).str.replace_all("-XX", ""))

                # Correct date format
                meta_corr = temp_meta.filter((pl.col(cd).str.len_chars() == 4) |
                                             (pl.col(cd).str.len_chars() == 7) |
                                             (pl.col(cd).str.len_chars() == 10))

                # INFO without incomplete dates ----------

                meta_compl = meta_corr.filter(pl.col(cd).str.len_chars() == len("YYYY-MM-DD"))

                meta_compl = meta_compl.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                meta_compl = meta_compl.filter((pl.col(cd) >= start) &
                                               (pl.col(cd) <= end))

                doc.write(f"\nLength of the metadata filtered on the basis of the collection dates"
                          f" WITHOUT incomplete collection dates considered: {len(meta_compl)}\n")

                # INFO with incomplete dates filled ----------

                meta_incompl = meta_corr.with_columns(pl.col(cd).map_elements(lambda x: x + "-01-01" if len(x) == 4 else x))
                meta_incompl = meta_incompl.with_columns(pl.col(cd).map_elements(lambda x: x + "-15" if len(x) == 7 else x))

                meta_incompl = meta_incompl.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                meta_incompl = meta_incompl.filter((pl.col(cd) >= start) &
                                                   (pl.col(cd) <= end))

                doc.write(f"\nLength of the metadata filtered on the basis of the collection dates"
                          f" WITH incomplete collection dates considered and FILLED WITH MISSING MONTH/DAY:"
                          f" {len(meta_incompl)}\n\n")

            elif start_date:

                # Dates formatting
                start = datetime.strptime(start_date, '%Y-%m-%d')

                temp_meta = meta.filter(pl.col(cd) != "Undefined")

                temp_meta = temp_meta.with_columns(pl.col(cd).str.replace_all("/", "-"))

                # Correct date format
                meta_corr = temp_meta.filter((pl.col(cd).str.len_chars() == 4) |
                                             (pl.col(cd).str.len_chars() == 7) |
                                             (pl.col(cd).str.len_chars() == 10))

                # INFO without incomplete dates ----------

                meta_compl = meta_corr.filter(pl.col(cd).str.len_chars() == len("YYYY-MM-DD"))

                meta_compl = meta_compl.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                meta_compl = meta_compl.filter(pl.col(cd) >= start)

                doc.write(f"\nLength of the metadata filtered on the basis of the collection dates"
                          f" WITHOUT incomplete collection dates considered: {len(meta_compl)}\n")

                # INFO with incomplete dates filled ----------

                meta_incompl = meta_corr.with_columns(
                    pl.col(cd).map_elements(lambda x: x + "-01-01" if len(x) == 4 else x))
                meta_incompl = meta_incompl.with_columns(
                    pl.col(cd).map_elements(lambda x: x + "-15" if len(x) == 7 else x))

                meta_incompl = meta_incompl.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                meta_incompl = meta_incompl.filter(pl.col(cd) >= start)

                doc.write(f"\nLength of the metadata filtered on the basis of the collection dates"
                          f" WITH incomplete collection dates considered and FILLED WITH MISSING MONTH/DAY:"
                          f" {len(meta_incompl)}\n\n")

            elif end_date:

                end = datetime.strptime(end_date, '%Y-%m-%d')

                temp_meta = meta.filter(pl.col(cd) != "Undefined")

                temp_meta = temp_meta.with_columns(pl.col(cd).str.replace_all("/", "-"))

                # Correct date format
                meta_corr = temp_meta.filter((pl.col(cd).str.len_chars() == 4) |
                                             (pl.col(cd).str.len_chars() == 7) |
                                             (pl.col(cd).str.len_chars() == 10))

                # INFO without incomplete dates ----------

                meta_compl = meta_corr.filter(pl.col(cd).str.len_chars() == len("YYYY-MM-DD"))

                meta_compl = meta_compl.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                meta_compl = meta_compl.filter(pl.col(cd) <= end)

                doc.write(f"\nLength of the metadata filtered on the basis of the collection dates"
                          f" WITHOUT incomplete collection dates considered: {len(meta_compl)}\n")

                # INFO with incomplete dates filled ----------

                meta_incompl = meta_corr.with_columns(
                    pl.col(cd).map_elements(lambda x: x + "-01-01" if len(x) == 4 else x))
                meta_incompl = meta_incompl.with_columns(
                    pl.col(cd).map_elements(lambda x: x + "-15" if len(x) == 7 else x))

                meta_incompl = meta_incompl.with_columns(pl.col(cd).str.to_date("%Y-%m-%d"))
                meta_incompl = meta_incompl.filter(pl.col(cd) <= end)

                doc.write(f"\nLength of the metadata filtered on the basis of the collection dates"
                          f" WITH incomplete collection dates considered and FILLED WITH MISSING MONTH/DAY:"
                          f" {len(meta_incompl)}\n\n")


# MAIN #################################################################################################################

parser = argparse.ArgumentParser(description="Metadata quality diagnostic.\n"
                                             "- Detection of human or any other required host sequences.\n"
                                             "- Detection of sequences collected in the required place.\n"
                                             "- Detection of sequences collected in the required dates.\n"
                                             "- Detection of sequences associated with incomplete collection dates."
                                             "\n"
                                             "- Detection of sequences associated with lacking information.\n",
                                 formatter_class=RawTextHelpFormatter)

# User metadata file
parser.add_argument("usr_metadata", type=str, help="Metadata file path.\nSequence ID, Collection date, Location "
                                                   "and Host columns are mandatory.\nIf you have multiple metadata files"
                                                   ", concatenate or merge them.\nMANDATORY FORMAT: tsv FILE.")

# Virus host
parser.add_argument("--host", type=str, help="Detection of sequences associated with the required host (optional).")

# Place
parser.add_argument("--place", type=str, help="Detection of sequences sampled in the required place (optional).\n")

# Start date
parser.add_argument("--start_date", type=str, help="Selection of sequences sampled after the required collection "
                                                   "date (optional).\nIf you want to crate a specific time span "
                                                   "add the --end_date optional parameter.")

# End date
parser.add_argument("--end_date", type=str, help="Selection of sequences sampled before the required collection "
                                                 "date (optional).\nIf you want to crate a specific time span "
                                                 "add the --start_date optional parameter.")

# Metadata columns
parser.add_argument("--meta_columns_list", type=str, help="If your metadata file doesn't match the GISAID and NCBI "
                                                          "formats, you MUST provide a list of columns as an extra input"
                                                          ". The list should follow this pattern:\n[\"sequence id "
                                                          "column name\", \"host column name\", \"place column name\","
                                                          " \"collection date column name\"].\n ATTENTION: if you have "
                                                          "missing columns, do maintain the order and use \"\" in its "
                                                          "place. For instance, if you lack the \"host\" column, your "
                                                          "list should look like this: [\"sequence id column name\", "
                                                          "\"\", \"place column name\", \"collection date column name"
                                                          "\"].\nATTENTION: the sequence ID column is mandatory.")

# Output folder path
parser.add_argument("--outfolder", type=str, help="Output folder path (optional).\nDefault out folder: current folder.")

args = parser.parse_args()

meta_cleaner_diagnostic(args.usr_metadata, args.outfolder, args.host, args.place, args.start_date, args.end_date,
                        args.meta_columns_list)
