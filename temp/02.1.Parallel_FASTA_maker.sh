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
# - w: white spaces filled with underscores through the file (y/n)
# ATTENTION: if you are working with NCBI FASTA files and you download
# not only the headers of the sequences but also associated
# descriptions/metadata, by setting this parameter as y the indexing
# process of your FASTA files could not work properly.
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
echo "- w: if you want the white spaces to be filled with underscores set as y, otherwise n."
echo "ATTENTION: if you are working with NCBI FASTA files and you download not only the headers of the"
echo "of the sequences but also associated descriptions/metadata, by setting this parameter as y the"
echo "indexing process of your FASTA files could not work properly."

# Command line inputs ----------
while getopts r:i:h:o:s:w: flag
do
    case "${flag}" in
        r) ref_folder_path=${OPTARG};;
        i) in_path=${OPTARG};;
        h) head_f_path=${OPTARG};;
        o) out_path=${OPTARG};;
        s) sif_folder=${OPTARG};;
        w) white_spaces=${OPTARG};;
    esac
done

echo "ref_folder_path = $ref_folder_path"
echo "in path = $in_path"
echo "headers folder = $head_f_path"
echo "out path = $out_path"
echo "sif folder = $sif_folder"
echo "filled white spaces = $white_spaces"

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

# FASTA files
if [[ $in_path == *.fasta ]]
then
        in_path=$in_path
        echo $in_path
elif [[ $in_path == */ ]]
then
        in_path=$in_path
else
        in_path=$in_path/
fi

temp_in_path=${ref_folder_path}${in_path}
echo $temp_in_path

# Unique FASTA file creation and old FASTAs zipping
if [[ $temp_in_path == */ ]]
then
	all_files=$(ls ${temp_in_path}*.fasta | wc -l)
	echo $all_files

	# If we have multiple FASTA files we create a unique one
	if [[ ${all_files} -gt 1 ]]
	then
		DIR=${temp_in_path}All_single_FASTAs
		if [ -d "$DIR" ]
		then
			temp_dir=$DIR
		else
			mkdir ${temp_in_path}All_single_FASTAs/
			temp_dir=${temp_in_path}All_single_FASTAs
		fi

		mv ${temp_in_path}*.fasta ${temp_dir}

		if [[ $white_spaces == "y" ]]
		then
			cat ${temp_dir}/*.fasta | sed 's/ /_/g' > ${temp_dir}/complete_FASTA.fasta

		elif [[ $white_spaces == "n" ]]
		then
			cat ${temp_dir}/*.fasta > ${temp_dir}/complete_FASTA.fasta
		fi

		mv ${temp_dir}/complete_FASTA.fasta ${temp_in_path}

		tar -czvf ${temp_dir}.tar.gz ${temp_dir}

		rm -r ${temp_dir}

		in_path=${in_path}complete_FASTA.fasta
		echo $in_path

	else
		in_path=${in_path}*.fasta
		echo $in_path
	fi
fi

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
pwd

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

# FASTA file white spaces removal and zip of the old FASTA file ----------
#ffile="$(basename -- $in_path .fasta)"
#
#fpath=${in_path%/*}
#echo "fasta file to parse: $ffile"
#
##echo "${in_path} sed 's/ /_/g' > ${fpath}/${ffile}_new.fasta"
#
##echo "zip ${in_path}.zip $in_path"
#
##echo "in_path=${fpath}/${ffile}_new.fasta"
##echo "fasta file to parse without white spaces: $in_path"
#
# faidx command line maker ----------
if [[ ${head_f_path} == *.txt ]]
then
	file="$(basename -- $head_f_path .txt)"
	echo "working on file: $file ----------"
	if test -f $out_path_dir${file}.fasta
	then
		echo "$out_path_dir${file}.fasta ALREADY EXISTS"
	else
		singularity exec --bind $ref_folder_path:$ref_folder_path ${sif_folder}samtools-1.17.sif samtools faidx ${in_path}${ffile} -r ${head_f_path} \"
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
			singularity exec --bind $ref_folder_path:$ref_folder_path ${sif_folder}samtools-1.17.sif samtools faidx $in_path${ffile} -r $file \
			2> $out_path_dir${filename}.log \
			1> ${out_path_dir}${filename}.fasta
		fi
	done
fi
