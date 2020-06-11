# $1 = reads
# $2 = transcriptome
# $3 = output_directory
# $4 = number_of_read_splitting_rounds
# $5 = make_transcript_bins [yes, no]
# $6 = delete_unnecessary_files [yes, no]
# $7 = reads_per_transcript_threshold
# $8 = scripts_folder

reads="$1"
ref="$2"
output_dir="$3"
n="$4"
make_transcript_bins="$5"
delete="$6"
threshold="$7"
scripts_folder="$8"

# variables which will be set by a python script

FA_tool="ram"
FA_PON="k10_w4"
FA_command='../ram/build/bin/ram -k 10 -w 4 $ref $reads > ${output_dir}/reads_selection/1st_split/alignment/$FA_output || ( echo >&2 && echo "ERROR in reads_selection.sh script: 1st_split ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1' #(will contain $reads $ref and $FA_output in it)

OA_tool="ram"
OA_PON="k10_w4"
OA_command='../ram/build/bin/ram -k 10 -w 4 ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_M ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_U > ${output_dir}/reads_selection/$current_split/alignment/${OA_output} || ( echo >&2 && echo "ERROR in reads_selection.sh script: ${current_split} ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1' #(will contain $OA_PRO_M $OA_PRO_U and $OA_output in it)


# check reads_file format

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


# check ref_file format ( we don't need its extension)

if echo $ref | grep -iq '\.fastq$' || echo $ref | grep -iq '\.fq$' ||  echo $ref | grep -iq '\.fasta$' || echo $ref | grep -iq '\.fa$' 
then

	if [ ! "$n" -ge 0 ]
	then
		echo >&2
		echo "ERROR in reads_selection.sh script: Incorrect argument: '$n'. Fourth argument (number_of_read_splitting_rounds) needs to be >= 0." >&2
		echo >&2
		exit 1
	fi

	if [ ! '(' "$make_transcript_bins" == 'yes' -o "$make_transcript_bins" == 'no' ')' ]
	then
		echo >&2
		echo "ERROR in reads_selection.sh script: Incorrect argument: '$make_transcript_bins'. Fifth argument (make_transcript_bins) needs to be 'yes' or 'no'." >&2
		echo >&2
		exit 1
	fi

	if [ ! '(' "$delete" == 'yes' -o "$delete" == 'no' ')' ]
	then
		echo >&2
		echo "ERROR in reads_selection.sh script: Incorrect argument: '$delete'. Sixth argument (delete_unnecessary_files) needs to be 'yes' or 'no'." >&2
		echo >&2
		exit 1
	fi

	if [ ! "$threshold" -ge 0 ]
	then
		echo >&2
		echo "ERROR in reads_selection.sh script: Incorrect argument: '$threshold'. Seventh argument (reads_per_transcript_threshold) needs to be >= 0." >&2
		echo >&2
		exit 1
	fi
	
	if [ ! '(' -e "$scripts_folder" -o -r "$scripts_folder" ')' ]
	then
		echo >&2
		echo "ERROR in reads_selection.sh script: Incorrect argument: '$scripts_folder'. Eight argument (scripts_folder) needs to accessible." >&2
		echo >&2
		exit 1
	fi

	# if number_of_read_splitting_rounds == 0 we don't have to do anything

	if [ "$n" -eq 0 ]
	then
		echo >&2
		echo "WARNING: number_of_read_splitting_rounds == 0. Exiting." >&2
		echo >&2
		exit 0
	fi



	if [ ! -e "$output_dir" ]
	then
		mkdir $output_dir || ( echo >&2 && echo "ERROR in reads_selection.sh script: Couldn't make output directory." >&2 && echo >&2 && exit 1 ) || exit 1

	else
		if [ -d "${output_dir}/reads_selection" ]
		then
			rm -r ${output_dir}/reads_selection
		fi

	fi


	# get reads_file_name and ref_file_name

	# remove path_to_file_folder (awk), then extension (sed)

	reads_name=$(echo $reads | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')
	ref_name=$(echo $ref | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')

	# make necessary folders

	mkdir ${output_dir}/reads_selection
	mkdir ${output_dir}/reads_selection/pickle_files


	# do the first split (if n >= 1)

	if [ "$n" -ge 1 ]
	then
		mkdir ${output_dir}/reads_selection/1st_split
		mkdir ${output_dir}/reads_selection/1st_split/alignment
		mkdir ${output_dir}/reads_selection/1st_split/split_files

		if [ -n "$FA_PON" ]
		then
			FA_output="${FA_tool}_alignment_${FA_PON}_${reads_name}_to_${ref_name}.paf"

		else
			FA_output="${FA_tool}_alignment_${reads_name}_to_${ref_name}.paf"
		fi


		#$FA_command

		#graphmap align --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r $ref -d $reads -o ${output_dir}/reads_selection/1st_split/alignment/$FA_output || ( echo >&2 && echo "ERROR in reads_selection.sh script: 1st_split ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1
		../ram/build/bin/ram -k 10 -w 4 $ref $reads > ${output_dir}/reads_selection/1st_split/alignment/$FA_output || ( echo >&2 && echo "ERROR in reads_selection.sh script: 1st_split ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1 #(will contain $reads $ref and $FA_output in it)

		#<FA_command> -r $ref -d $reads -o ${output_dir}/reads_selection/1st_split/alignment/$FA_output || ( echo >&2 && echo "ERROR in reads_selection.sh script: 1st_split ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1

		python ${scripts_folder}/split_mapped_and_unmapped.py ${output_dir}/reads_selection/1st_split/alignment/$FA_output $reads ${output_dir}/reads_selection/1st_split/split_files ${output_dir}/reads_selection/pickle_files || ( echo >&2 && echo "ERROR in reads_selection.sh script: 1st SPLITTING failed!" >&2 && echo >&2 && exit 1 ) || exit 1

		OA_PRO_M="${reads_name}_1st_split_only_mapped${extension}"
		OA_PRO_U="${reads_name}_1st_split_only_unmapped${extension}"

		if [ "$delete" == 'yes' ]
		then
			rm -r ${output_dir}/reads_selection/1st_split/alignment
		fi

	fi


	if [ "$n" -ge 2 ]
	then
		split_number=2

		while [ "$split_number" -le "$n" ]
		do
			if [ "$split_number" -eq 2 ]
			then
				current_split="2nd_split"
				previous_split="1st_split"

			elif [ "$split_number" -eq 3 ]
			then
				current_split="3rd_split"
				previous_split="2nd_split"

			elif [ "$split_number" -eq 4 ]
			then
				current_split="4th_split"
				previous_split="3rd_split"

			else
				current_split="${split_number}th_split"
				previous_split="$((split_number - 1))th_split"
			fi


			mkdir ${output_dir}/reads_selection/$current_split
			mkdir ${output_dir}/reads_selection/$current_split/alignment
			mkdir ${output_dir}/reads_selection/$current_split/split_files


			if [ -n "$OA_PON" ]
			then
				OA_output="${OA_tool}_alignment_${OA_PON}_${previous_split}_unmapped_to_mapped.paf"

			else
				OA_output="${OA_tool}_alignment_${previous_split}_unmapped_to_mapped.paf"
			fi


			#$OA_command

			#graphmap align --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_M -d ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_U -o ${output_dir}/reads_selection/$current_split/alignment/${OA_output} || ( echo >&2 && echo "ERROR in reads_selection.sh script: ${current_split} ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1
			../ram/build/bin/ram -k 10 -w 4 ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_M ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_U > ${output_dir}/reads_selection/$current_split/alignment/${OA_output} || ( echo >&2 && echo "ERROR in reads_selection.sh script: ${current_split} ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1 #(will contain $OA_PRO_M $OA_PRO_U and $OA_output in it)

			#<OA_command> -r ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_M -d ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_U -o ${output_dir}/reads_selection/$current_split/alignment/${OA_output} || ( echo >&2 && echo "ERROR in reads_selection.sh script: ${current_split} ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1

			python ${scripts_folder}/split_mapped_and_unmapped.py ${output_dir}/reads_selection/${current_split}/alignment/$OA_output ${output_dir}/reads_selection/${previous_split}/split_files/$OA_PRO_U ${output_dir}/reads_selection/${current_split}/split_files ${output_dir}/reads_selection/pickle_files || ( echo >&2 && echo "ERROR in reads_selection.sh script: ${current_split} SPLITTING failed!" >&2 && echo >&2 && exit 1 ) || exit 1


			OA_PRO_M="${reads_name}_${current_split}_only_mapped${extension}"
			OA_PRO_U="${reads_name}_${current_split}_only_unmapped${extension}"


			if [ "$delete" == 'yes' ]
			then
				rm -r ${output_dir}/reads_selection/${current_split}/alignment

				if [ "$make_transcript_bins" == 'yes' ]
				then
					rm -r ${output_dir}/reads_selection/${previous_split}
				fi
			fi

			split_number=$((split_number + 1))

		done
	fi




	if [ "$make_transcript_bins" == 'yes' ]
	then
		mkdir ${output_dir}/transcript_bins

        #removed --delete_pickle_files!!!
		python ${scripts_folder}/create_transcript_bins.py ${output_dir}/reads_selection/pickle_files $reads ${output_dir}/transcript_bins --threshold $threshold || ( echo >&2 && echo "ERROR in reads_selection.sh script: CREATING TRANSCRIPT BINS FAILED!" >&2 && echo >&2 && exit 1 ) || exit 1

	elif [ "$make_transcript_bins" == 'no' ]
	then
		mkdir ${output_dir}/merged_files

		for file in ${output_dir}/reads_selection/*_split/split_files/*_only_mapped*
		do
			cat $file
		done > ${output_dir}/merged_files/"${reads_name}${extension}"
	fi

	if [ "$delete" == 'yes' ]
	then
		rm -r ${output_dir}/reads_selection
	fi

else
	echo "ERROR in reads_selection.sh script: Incorrect argument: '$ref'. Second argument (transcriptome_file) needs to be FASTA or FASTQ file." >&2 | exit 1
fi




