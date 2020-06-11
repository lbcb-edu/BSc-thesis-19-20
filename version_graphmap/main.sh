reads=$1
ref=$2
output_dir=$3
bins=$4
n=3

if test "$#" -ne 4
then
        echo "Illegal number of arguments. This script is run with path to reads, path to the  reference, output directory and yes/no depending on if transcript bins should be created"
        exit 1
fi

ref_name=$(echo $ref | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')
reads_name=$(echo $reads | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')

if echo $reads | grep -iq '\.fastq$'
then
        extension='.fastq'

elif echo $reads | grep -iq '\.fq$'
then
        extension='.fq'

elif echo $reads | grep -iq '\.fasta$'
then
        extension='.fasta'

elif echo $reads | grep -iq '\.fa$'
then
        extension='.fa'

else
        echo >&2
        echo "ERROR in reads_selection.sh script: Incorrect argument: '$reads'. First argument (reads_file) needs to be FASTA or FASTQ file." >&2
        echo >&2
        exit 1
fi

if $bins -eq "yes"
then 
        ./reads_selection.sh $reads $ref $output_dir $n yes no 3 .
        polish_reads=${output_dir}/transcript_bins/"${ref_name}${extension}"
else
        ./reads_selection.sh $reads $ref $output_dir $n no no 3 .
        polish_reads=${output_dir}/merged_files/"${reads_name}${extension}"
fi

./polish.sh $polish_reads $ref $output_dir $n no yes


P_output="polished_${n}x_${reads_name}.fasta"
PR_FN=$(echo $P_output | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')
FINAL_tool="graphmap"
FINAL_PON="k6_A0.5"

FINAL_output="${FINAL_tool}_alignment_${FINAL_PON}_${PR_FN}_to_${ref_name}.sam"

python clustering.py ${output_dir}/final_alignment/${FINAL_output} $ref



