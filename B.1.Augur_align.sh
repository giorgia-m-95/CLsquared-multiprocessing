# INFO #################################################################
#
# INPUTS ---------------------------------------------------------------
# - i: FASTA files to align path
# - r: Reference sequence path
# - o: Out folder path
# - t: Number of threads
# ----------------------------------------------------------------------
#
# DESCRIPTION ----------------------------------------------------------
# Automated augur align command line for multiple FASTA files to align
# versus the pathogen reference sequence.
# ----------------------------------------------------------------------
#
# OUTPUTS --------------------------------------------------------------
# - FASTA files with starting sequences aligned versus the pathogen
#   reference sequence with associated .log files.
# ----------------------------------------------------------------------
#
# ######################################################################

echo "- i: path of the FASTA files that you want to align"
echo "- r: pathogen reference sequence path"
echo "- o: out folder path"
echo "- t: number of threads"

# Command line inputs ----------
while getopts i:r:o:t: flag
do
    case "${flag}" in
        i) in_path=${OPTARG};;
        r) ref_path=${OPTARG};;
        o) out_path=${OPTARG};;
        t) n_threads=${OPTARG};;
    esac
done


# Inputs format control ----------
if [[ $in_path == *.fasta ]]
then
	in_path=$in_path
elif [[ $in_path == */ ]]
then
	in_path=$in_path
else
	in_path=$in_path/
fi

if [[ $out_path == */ ]]
then
	out_path=$out_path
else
	out_path=$out_path/
fi


DIR=${out_path}CL_squared_step_alignment/
if [ -d "$DIR" ]
then
	out_path_dir=$DIR
else
	mkdir ${out_path}CL_squared_alignment/
	out_path_dir=${out_path}CL_squared_alignment/
fi


# augur align command line maker ----------
if [[ $in_path == *.fasta ]]
then
	file="$(basename -- $in_path .fasta)"
	if test -f $out_path_dir${file}_aligned.fasta.log
	then
		echo $out_path_dir${file}_aligned.fasta.log ALREADY EXISTS
	else
		out_file_name="$(basename -- $in_path .fasta)"
		augur align --sequences $in_path --reference-sequence $ref_path --output $out_path_dir${out_file_name}_aligned.fasta --nthreads $n_threads
	fi
else
	for file in $in_path*.fasta
	do
		filename="$(basename -- $file .fasta)"
		if test -f $out_path_dir${filename}_aligned.fasta.log
		then
			echo $out_path_dir${filename}_aligned.fasta.log ALREADY EXISTS
		else
			out_file_name="$(basename -- $file .fasta)"
			augur align --sequences $file --reference-sequence $ref_path --output $out_path_dir${out_file_name}_aligned.fasta --nthreads $n_threads
		fi
	done
fi
