# INFO #################################################################
#
# INPUTS ---------------------------------------------------------------
# - r: folder to bind. All the other paths will refer to this folder.
# - i: path of the FASTA file to split in multiple files (with respect
#      to r)
# - h: path of the folder containing the .txt files with the groups of
#      headers (with respect to r)
# - o: out folder path (with respect to r)
# - s: path of the faidx sif folder to bind (with respect to r)
# ----------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------
# The script splits a huge FASTA file into multiple files on the basis
# of the .txt files containing the desidered groups of headers. The
# resulting FASTA files will be used to parallelize the CLsquared
# pipeline steps.
# The faidx tool is exploited.
# It is necessary to indicate the absolute path of the folder
# containing the container to bind.
# ----------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------
# - FASTA files obtained from the split of the big starting FASTA.
# ----------------------------------------------------------------------
#
# ######################################################################

echo "- r: folder to bind. All the other paths will refer to this folder."
echo "- i: path of the FASTA file that you want to parse (with respect to r)"
echo "- h: path of the folder containing the .txt files with the groups of headers (with respect to r)"
echo "- o: out folder path (with respect to r)"
echo "- s: path of the faidx sif folder to bind (with respect to r)"

# Command line inputs ----------
while getopts r:i:h:o:s: flag
do
    case "${flag}" in
		r) ref_folder_path=${OPTARG};;
        i) in_path=${OPTARG};;
        h) head_f_path=${OPTARG};;
        o) out_path=${OPTARG};;
		s) sif_folder=${OPTARG};;
    esac
done

echo "ref_folder_path = $ref_folder_path"
echo "in path = $in_path"
echo "headers folder = $head_f_path"
echo "out path = $out_path"
echo "sif folder = $sif_folder"

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
echo "Folder to bind: $ref_folder_path"

# Headers files
if [[ $head_f_path == *.txt ]]
then
	head_f_path=$head_f_path
elif [[ $head_f_path == */ ]]
then
	head_f_path=$head_f_path
else
	head_f_path=$head_f_path/
fi
echo "Headers folder: $head_f_path"

# Out path
if [[ $out_path == */ ]]
then
	out_path=$out_path
else
	out_path=$out_path/
fi
echo "Out folder: $out_path"

# sif folder
if [[ $sif_folder == */ ]]
then
	sif_folder=$sif_folder
else
	sif_folder=$sif_folder/
fi
echo "Sif folder: $sif_folder"

cd ${ref_folder_path}

# Out folder maker ----------
DIR=${out_path}CL_squared_parallel_FASTAs/
if [ -d "$DIR" ]
then
	 out_path_dir=$DIR
 else
	 mkdir ${out_path}CL_squared_parallel_FASTAs/
	 out_path_dir=${out_path}CL_squared_parallel_FASTAs/
fi
echo "Out path directory: $out_path_dir"

## FASTA file white spaces removal and zip of the old FASTA file ----------
#ffile="$(basename -- $in_path .fasta)"
#fpath=${in_path%/*}
#echo "fasta file to parse: $ffile"

#echo "${in_path} sed 's/ /_/g' > ${fpath}/${ffile}_new.fasta"

#echo "zip ${in_path}.zip $in_path"

#echo "in_path=${fpath}/${ffile}_new.fasta"
#echo "fasta file to parse without white spaces: $in_path"

# faidx command line maker ----------
if [[ ${head_f_path} == *.txt ]]
then
	file="$(basename -- $head_f_path .txt)"
	echo "working on file: $file ----------"
	if test -f $out_path_dir${file}.fasta
	then
		echo "$out_path_dir${file}.fasta ALREADY EXISTS"
	else
		singularity exec --bind $ref_folder_path:$ref_folder_path ${sif_folder}samtools-1.17.sif samtools faidx $in_path -r $head_f_path \"
		2> $out_path_dir${file}.log \
		1> ${out_path_dir}${file}.fasta
	fi
else
	for file in ${head_f_path}*.txt
	do
		filename="$(basename -- ${file} .txt)"
		echo "working on file: ${filename}  ----------"
		if test -f $out_path_dir${filename}.fasta
		then
			echo "$out_path_dir${filename}.fasta ALREADY EXISTS"
		else
			singularity exec --bind $ref_folder_path:$ref_folder_path ${sif_folder}samtools-1.17.sif samtools faidx $in_path -r $file \
			2> $out_path_dir${filename}.log \
			1> ${out_path_dir}${filename}.fasta
		fi
	done
fi
