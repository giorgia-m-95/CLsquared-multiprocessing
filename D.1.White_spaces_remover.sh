# INFO #################################################################
#
# INPUTS ---------------------------------------------------------------
# - r: folder to bind. All the other paths will refer to this folder.
# - i: path of the FASTA file to split in multiple files (with respect
#      to r)
# - o: out folder path (with respect to r)
# ----------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------
# The script substitutes the white spaces in the sequence headers of the
# FASTA files with underscores.
# Please pay attention when working with NCBI sequences.
# ----------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------
# - FASTA files with white spaces substituted with underscores
# ----------------------------------------------------------------------
#
# ######################################################################

echo "- r: folder to bind. All the other paths will refer to this folder."
echo "- i: path of the FASTA file that you want to parse (with respect to r)"
echo "- o: out folder path (with respect to r)"

# Command line inputs ----------
while getopts r:i:o: flag
do
    case "${flag}" in
        r) ref_folder_path=${OPTARG};;
        i) in_path=${OPTARG};;
        o) out_path=${OPTARG};;
    esac
done

echo "ref_folder_path = $ref_folder_path"
echo "in path = $in_path"
echo "out path = $out_path"

# Inputs format control ----------
# Refrence folder
if [[ $ref_folder_path == */ ]]
then
        ref_folder_path=$ref_folder_path
else
        ref_folder_path=$ref_folder_path/
fi

if [[ $ref_folder_path == /* ]]
then
        ref_folder_path=$ref_folder_path
else
        ref_folder_path=/${ref_folder_path}
fi

# FASTA files
if [[ $in_path == *.fasta ]]
then
        in_path=$in_path
        echo $in_path
elif [[ $in_path == */ ]]
then
        in_path=$in_path*.fasta
else
        in_path=$in_path/*.fasta
fi

temp_in_path=${ref_folder_path}${in_path}

# Out path
if [[ $out_path == */ ]]
then
        out_path=$out_path
else
        out_path=$out_path/
fi

# Out folder maker ----------
DIR=${ref_folder_path}${out_path}FASTAs_no_white_spaces/

if [ -d "$DIR" ]
then
         out_path_dir=$DIR
else
         mkdir ${ref_folder_path}${out_path}FASTAs_no_white_spaces/
         out_path_dir=${ref_folder_path}${out_path}FASTAs_no_white_spaces/
fi

echo ${out_path_dir}

for f in ${temp_in_path}
do
	name=$(basename ${f} .fasta)
        echo ${name}
        cat ${f} | sed 's/ /_/g' > ${out_path_dir}${name}_no_spaces_FASTA.fasta
done
