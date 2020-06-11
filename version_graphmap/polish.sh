# $1 = reads
# $2 = reference
# $3 = output_directory
# $4 = number_of_polishing_rounds
# $5 = delete_unnecessary_files [yes, no]
# $6 = do_final_alignment [yes, no]

reads="$1"
ref="$2"
output_dir="$3"
n="$4"
delete="$5"
do_final_alignment="$6"

# variables which will be set by a python script

FINAL_tool="graphmap"
FINAL_PON="k6_A0.5"
FINAL_command='graphmap align --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r $ref -d ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$polished_reads -o ${output_dir}/final_alignment/$FINAL_output || ( echo >&2 && echo "ERROR in polish.sh script: FINAL ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1' #(will contain $polished_reads $ref and $FINAL_output in it)

OV_tool="graphmap"
OV_PON="k6_A0.5"

OV_command_1='graphmap align -x overlap --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r $reads_to_polish -d $reads_to_polish -o ${output_dir}/polishing/1x/base_to_base_overlap/$OV_output || ( echo >&2 && echo "ERROR in polish.sh script: 1st OVERLAP failed!" >&2 && echo >&2 && exit 1 ) || exit 1' #(will contain $reads_to_polish and $OV_output in it)

OV_command_2='graphmap align -x overlap --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish -d ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish -o ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap/$OV_output || ( echo >&2 && echo "ERROR in polish.sh script: $current_round OVERLAP failed!" >&2 && echo >&2 && exit 1 ) || exit 1'

P_tool="racon" #(doesn't get used)
P_PON="" #(doesn't get used)

P_command_1='racon -u -q 0 -e 0.7 -m 4 -x -4 -g -4 -t 4 -f $reads_to_polish ${output_dir}/polishing/1x/base_to_base_overlap/$OV_output $reads_to_polish > ${output_dir}/polishing/1x/polished_files/$P_output || ( echo >&2 && echo "ERROR in polish.sh script: 1st POLISH failed!" >&2 && echo >&2 && exit 1 ) || exit 1' #(will contain $reads_to_polish, $OV_output and $P_output in it)

P_command_2='racon -u -q 0 -e 0.7 -m 4 -x -4 -g -4 -t 4 -f ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap/$OV_output ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish > ${output_dir}/polishing/${polishing_round}x/polished_files/$P_output || ( echo >&2 && echo "ERROR in polish.sh script: $current_round POLISH failed!" >&2 && echo >&2 && exit 1 ) || exit 1'

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
	echo "ERROR in polish.sh script: Incorrect argument: '$reads'. First argument (reads_file) needs to be FASTA or FASTQ file." >&2
	echo >&2
	exit 1
fi

# check ref_file format ( we don't need its extension)

if echo $ref | grep -iq '\.fastq$' || echo $ref | grep -iq '\.fq$' ||  echo $ref | grep -iq '\.fasta$' || echo $ref | grep -iq '\.fa$' 
then

	if [ ! "$n" -ge 0 ]
	then
		echo >&2
		echo "ERROR in polish.sh script: Incorrect argument: '$n'. Fourth argument (number_of_polishing_rounds) needs to be >= 0." >&2
		echo >&2
		exit 1
	fi

	if [ ! '(' "$delete" == 'yes' -o "$delete" == 'no' ')' ]
	then
		echo >&2
		echo "ERROR in polish.sh script: Incorrect argument: '$delete'. Fifth argument (delete_unnecessary_files) needs to be 'yes' or 'no'." >&2
		echo >&2
		exit 1
	fi

	if [ ! '(' "$do_final_alignment" == 'yes' -o "$do_final_alignment" == 'no' ')' ]
	then
		echo >&2
		echo "ERROR in polish.sh script: Incorrect argument: '$do_final_alignment'. Sixth argument (do_final_alignment) needs to be 'yes' or 'no'." >&2
		echo >&2
		exit 1
	fi

	# if number_of_polishing_rounds == 0 we don't have to do anything

	if [ "$n" -eq 0 ]
	then
		echo >&2
		echo "WARNING in polish.sh script: 'number_of_polishing_rounds' == 0. Exiting." >&2
		echo >&2
		exit 0
	fi



	if [ ! -e "$output_dir" ]
	then
		mkdir $output_dir || ( echo >&2 && echo "ERROR in polish.sh script: Couldn't make output directory." >&2 && echo >&2 && exit 1 ) || exit 1

	else
		if [ -d "${output_dir}/polishing" ]
		then
			rm -r ${output_dir}/polishing
		fi

		if [ -d "${output_dir}/final_alignment" ]
		then
			rm -r ${output_dir}/final_alignment
		fi
	fi

	# get reads_file_name and ref_file_name

	# remove path_to_file_folder (awk), then extension (sed)

	reads_name=$(echo $reads | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')
	ref_name=$(echo $ref | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')

	# make necessary folders

	mkdir ${output_dir}/polishing


	# do the first polishing round (if n >= 1)

	if [ "$n" -ge 1 ]
	then
		mkdir ${output_dir}/polishing/1x
		mkdir ${output_dir}/polishing/1x/base_to_base_overlap
		mkdir ${output_dir}/polishing/1x/polished_files

		reads_to_polish=$reads

		# RTP_FN = reads_to_polish_file_name
		RTP_FN=$reads_name

		if [ -n "$OV_PON" ]
		then
			OV_output="overlap_${OV_tool}_${OV_PON}_${RTP_FN}_to_${RTP_FN}.sam"

		else
			OV_output="overlap_${OV_tool}_${RTP_FN}_to_${RTP_FN}.sam"
		fi

		#$OV_command_1

		graphmap align -x overlap --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r $reads_to_polish -d $reads_to_polish -o ${output_dir}/polishing/1x/base_to_base_overlap/$OV_output || ( echo >&2 && echo "ERROR in polish.sh script: 1st OVERLAP failed!" >&2 && echo >&2 && exit 1 ) || exit 1

		#<OV_command> -r $reads_to_polish -d $reads_to_polish -o ${output_dir}/polishing/1x/base_to_base_overlap/$OV_output || ( echo >&2 && echo "ERROR in polish.sh script: 1st OVERLAP failed!" >&2 && echo >&2 && exit 1 ) || exit 1


		# uncomment this section if you want to have P_tool and P_PON in P_output

		#if [ -n "$P_PON" ]
		#then
		#	P_output="polished_1x_${P_tool}_${P_PON}_${reads_name}.fasta"
		#
		#else
		#	P_output="polished_1x_${P_tool}_${reads_name}.fasta"
		#fi

		P_output="polished_1x_${reads_name}.fasta"

		#$P_command_1

		racon -u -q 0 -e 0.7 -m 4 -x -4 -g -4 -t 4 -f $reads_to_polish ${output_dir}/polishing/1x/base_to_base_overlap/$OV_output $reads_to_polish > ${output_dir}/polishing/1x/polished_files/$P_output || ( echo >&2 && echo "ERROR in polish.sh script: 1st POLISH failed!" >&2 && echo >&2 && exit 1 ) || exit 1

		#<P_command> -f $reads_to_polish ${output_dir}/polishing/1x/base_to_base_overlap/$OV_output $reads_to_polish > ${output_dir}/polishing/1x/polished_files/$P_output || ( echo >&2 && echo "ERROR in polish.sh script: 1st POLISH failed!" >&2 && echo >&2 && exit 1 ) || exit 1

		if [ "$delete" == 'yes' ]
		then
			rm -r ${output_dir}/polishing/1x/base_to_base_overlap
		fi

	fi


	if [ "$n" -ge 2 ]
	then
		polishing_round=2

		while [ "$polishing_round" -le "$n" ]
		do
			if [ "$polishing_round" -eq 2 ]
			then
				current_round="2nd"

			elif [ "$polishing_round" -eq 3 ]
			then
				current_round="3rd"

			else
				current_round="${polishing_round}th"
			fi


			mkdir ${output_dir}/polishing/${polishing_round}x
			mkdir ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap
			mkdir ${output_dir}/polishing/${polishing_round}x/polished_files


			reads_to_polish=$P_output

			# reads_to_polish_file_name
			RTP_FN=$(echo $P_output | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')


			if [ -n "$OV_PON" ]
			then
				OV_output="overlap_${OV_tool}_${OV_PON}_${RTP_FN}_to_${RTP_FN}.sam"

			else
				OV_output="overlap_${OV_tool}_${RTP_FN}_to_${RTP_FN}.sam"
			fi

			#$OV_command_2

			graphmap align -x overlap --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish -d ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish -o ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap/$OV_output || ( echo >&2 && echo "ERROR in polish.sh script: $current_round OVERLAP failed!" >&2 && echo >&2 && exit 1 ) || exit 1

			#<OV_command> -r ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish -d ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish -o ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap/$OV_output || ( echo >&2 && echo "ERROR in polish.sh script: $current_round OVERLAP failed!" >&2 && echo >&2 && exit 1 ) || exit 1

			# uncomment this section if you want to have P_tool and P_PON in P_output

			#if [ -n "$P_PON" ]
			#then
			#	P_output="polished_${polishing_round}x_${P_tool}_${P_PON}_${reads_name}.fasta"
			#
			#else
			#	P_output="polished_${polishing_round}x_${P_tool}_${reads_name}.fasta"
			#fi

			P_output="polished_${polishing_round}x_${reads_name}.fasta"

			#$P_command_2

			racon -u -q 0 -e 0.7 -m 4 -x -4 -g -4 -t 4 -f ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap/$OV_output ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish > ${output_dir}/polishing/${polishing_round}x/polished_files/$P_output || ( echo >&2 && echo "ERROR in polish.sh script: $current_round POLISH failed!" >&2 && echo >&2 && exit 1 ) || exit 1

			#<P_command> -f ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap/$OV_output ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$reads_to_polish > ${output_dir}/polishing/${polishing_round}x/polished_files/$P_output || ( echo >&2 && echo "ERROR in polish.sh script: $current_round POLISH failed!" >&2 && echo >&2 && exit 1 ) || exit 1



			if [ "$delete" == 'yes' ]
			then
				rm -r ${output_dir}/polishing/${polishing_round}x/$((polishing_round - 1))x_to_$((polishing_round - 1))x_overlap
				rm -r ${output_dir}/polishing/$((polishing_round - 1))x
			fi

			polishing_round=$((polishing_round + 1))

		done
	fi



	if [ "$do_final_alignment" == 'yes' ]
	then
		mkdir ${output_dir}/final_alignment

		# reads_to_align
		polished_reads=$P_output

		# polished_reads_file_name
		PR_FN=$(echo $P_output | awk -F / '{print $(NF)}' | sed -r 's/(.*)\..*/\1/')

		if [ -n "$FINAL_PON" ]
		then
			FINAL_output="${FINAL_tool}_alignment_${FINAL_PON}_${PR_FN}_to_${ref_name}.sam"

		else
			FINAL_output="${FINAL_tool}_alignment_${PR_FN}_to_${ref_name}.sam"
		fi

		#$FINAL_command

		graphmap align --freq-percentile 1.0 --double-index --minimizer-window 1 -k 6 --ambiguity 0.5 --auto-rebuild-index -r $ref -d ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$polished_reads -o ${output_dir}/final_alignment/$FINAL_output || ( echo >&2 && echo "ERROR in polish.sh script: FINAL ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1

		#<FINAL_command> -r $ref -d ${output_dir}/polishing/$((polishing_round - 1))x/polished_files/$polished_reads -o ${output_dir}/final_alignment/$FINAL_output || ( echo >&2 && echo "ERROR in polish.sh script: FINAL ALIGNMENT failed!" >&2 && echo >&2 && exit 1 ) || exit 1
	fi


else
	echo >&2
	echo "ERROR in polish.sh script: Incorrect argument: '$ref'. Second argument (reference file) needs to be FASTA or FASTQ file." >&2
	echo >&2
	exit 1
fi




