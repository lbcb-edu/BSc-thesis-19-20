# To find out what this script does run the script with '-h' or '--help'.

import re
import os
import pickle
import argparse
import errno

# Functions that are used to check input formats and make output directories if needed:


def output_dir(x):
    if not os.path.exists(x):
        try:
            os.makedirs(x, 0o700)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise

    return x


def pickle_dir(x):
    if not os.path.exists(x):
        try:
            os.makedirs(x, 0o700)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise

    return x


def paf_file(x):
    if not x.lower().endswith('.paf'):
        raise argparse.ArgumentTypeError("'paf_file' argument needs to be PAF file.")

    return x


def fa_or_fq_file(x):
    # we use extension variable to store if the reads_to_split file is FASTA or FASTQ
    global extension

    fastq_extensions = ('.fastq', '.fq')
    fasta_extensions = ('.fasta', '.fa')

    if x.lower().endswith(fastq_extensions):
        extension = '.fastq'

    elif x.lower().endswith(fasta_extensions):
        extension = '.fasta'

    else:
        raise argparse.ArgumentTypeError("'reads_to_split' argument needs to be FASTA or FASTQ file.")

    return x


# define arguments

# Argparse is a module that allows us to parse command line arguments easier:
# -> quick intro: https://docs.python.org/3/howto/argparse.html

parser = argparse.ArgumentParser(
    description=
    """
    This script takes as input: reads which were mapped to transcript/already mapped reads and PAF alignment file. 
    Using both it splits unmapped from mapped reads. Since there will be multiple splits, the script will track to 
    witch transcript do reads map by creating pickle files for each split. The files will be stored in separate 
    pickle directory.
    
    Pickle files will contain dictionaries of following type:
        ### transcript_bins.p:
            transcript -> list of reads that mapped to it
        ### read_mappings.p:
            read -> list of transcripts it mapped to.
                
    All dictionaries will be merged together and used to split reads into transcript bins with
    create_transcript_bins.py script.
    
    Naming conventions for output folders/files: 
        ### output folder names: 
            1st_split, 2nd_split, ... , <n>th_split
        ### output file names:
            original_full_reads_file_name + <split_number>th_split + only_mapped/only_unmapped + .fastq/.fasta.
    """
)

parser.add_argument(
    'paf_file', type=paf_file,

    help=
    """
    Alignment file in PAF format.
    
    NOTE:
        Script will check if the file ends with '.paf' or any of its variations.
        If not 'ArgumentTypeError' will be raised.
    """
)

parser.add_argument(
    'reads_to_split', type=fa_or_fq_file,

    help=
    """
    Reads to split in FASTQ/FASTA format.
    
    NOTE: 
        Script will check if the file ends with '.fastq', '.fasta', '.fq', '.fa' or any of its variations.
        If not 'ArgumentTypeError' will be raised.
    
        
    Since reads_to_split input file is either full reads' original file or one of the already split reads'
    only unmapped file, the script keeps track of split number by looking for 1st_split, 2nd_split, ... , <n>th_split
    in 'reads_to_split' file name or assumes it's the 1st split if it can't find anything.
    
    Example: 
        If the file name contains 1st_split then it means we are doing the 2nd split.
        If it can't find anything then it means it's the 1st_split.
    """
)

parser.add_argument(
    'output_directory', type=output_dir,

    help=
    """
    Output directory path.
    
    NOTE: 
        This directory isn't the top directory of output files. The script will just write mapped and unmapped read
        files in it. It should be the split's folder.
    """
)

parser.add_argument(
    'pickle_directory', type=pickle_dir,

    help=
    """
    Pickle files directory path.
    
    NOTE:
        This is just the top directory of pickle files. If it doesn't exist the script will try to create it.
        The script will find out which split it is from the reads_to_split file name, create current split's folder
        in pickles directory and place pickle files in it.
    """
)

args = parser.parse_args()

# get the name of the file (without .fastq or .fq / .fasta or .fa at the end)
parts = re.split('/', args.reads_to_split)
parts2 = re.split('\.', parts[-1])

file_name = parts2[0]

# depending on the name of reads_to_split file create current split's output folder and output files' paths

# which_split variable tells us which split it is according to the reads_to_split file name
# split_number contains previous split's number as its first element (it is a list of strings)

# program looks for 1st_split, 2nd_split, 3rd_split or <number>th_split in the name of reads_to_split file
# example: if the file name contains 1st_split then it means we are doing the 2nd split
# if it can't find anything then it means it's the 1st_split

which_split = ''
split_number = ''

# first we figure out which split it is and edit file_name to match the current split's name
if re.findall('1st_split', file_name):

    file_name = file_name.replace('1st_split_only_unmapped', '2nd_split')
    which_split = '2nd_split'

elif re.findall('2nd_split', file_name):

    file_name = file_name.replace('2nd_split_only_unmapped', '3rd_split')
    which_split = '3rd_split'

elif re.findall('3rd_split', file_name):

    file_name = file_name.replace('3rd_split_only_unmapped', '4th_split')
    which_split = '4th_split'

elif re.findall('[0-9]+th_split', file_name):

    tmp = re.findall('[0-9]+th_split', file_name)
    split_number = re.split('th', tmp[0])

    file_name = file_name.replace(tmp[0] + '_only_unmapped', str(int(split_number[0]) + 1) + 'th_split')
    which_split = str(int(split_number[0]) + 1) + 'th_split'

else:
    file_name = file_name + '_1st_split'
    which_split = '1st_split'


# create current split's output directory

if not os.path.exists(args.output_directory + '/' + which_split):
    try:
        os.makedirs(args.output_directory + '/' + which_split, 0o700)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

# set output files' paths
out_path_only_mapped = args.output_directory + '/' + file_name + '_only_mapped' + extension
out_path_only_unmapped = args.output_directory + '/' + file_name + '_only_unmapped' + extension


# read paf file

file_paf = open(args.paf_file, 'r', encoding='utf-8')


# Dictionaries:
#   -> current_split_dict_paf: (dictionary of type: read -> list of transcripts to which it mapped)
#   -> previous_split_dict_paf (dictionary of pafe type but made in previous split)
#   -> mapping_to_transcripts (dictionary of type: transcript -> list of reads which mapped to it)

#   -> comprehensive explanation is given in the next comment

current_split_dict_paf = {}
previous_split_dict_paf = {}
mapping_to_transcripts = {}

line = file_paf.readline()

# Next part does the following:
# -> For the current split creates a dictionary: transcript -> list of reads that mapped to it.
#    The dictionary is stored in <split_number><st/nd/rd/th>_split/transcript_bins.p file in pickle directory.
#    All split's dictionaries of this type will be merged into 1 dictionary with a separate script.
# -> For the current split creates a dictionary: read -> list of transcripts to which it mapped.
#    This dictionary will be stored under <split_number><st/nd/rd/th>_split/read_mappings.p file in pickle directory.
#    If it's the 1st split just create the dictionary.
#    For every succeeding split get the previous split's dictionary and use it to form the current dictionary.

# -> This is necessary due to the following:
#    After the first mapping, we are mapping unmapped reads to mapped ones. On account of that,
#    references in .paf file will be id's of previously mapped reads. To determine to which transcript bin
#    our read belongs ( or to witch transcript it 'maps' to ), we will use the previous split's dictionary.
#    From it, we can get the transcripts to which previously mapped read maps. Then we can say our read maps
#    to the pafe transcripts.

if which_split == '1st_split':
    while line != '':
        parts = re.split('\s+', line)

        # parts[0] = id of the mapped read

        # parts[1] = number which tells us the following:
        #            0 - read mapped
        #            4 - read didn't map
        #            16 - read mapped to a reverse complement

        # parts[2] = name of the reference the read mapped to

        read_id = parts[0]
        name_of_the_reference = parts[5]

        # We will read each line of PAF file and if the read mapped, add (key, value) pair to both
        # 'mapping_to_transcripts' and 'current_split_dict_paf'. In 'mapping_to_transcripts' keys are transcript names
        # and values lists of read_ids that mapped to the transcript. In 'current_split_dict_paf' keys are read_ids and
        # values lists of transcripts they mapped to. We have two situations:
        # 1) first appearance of a key:
        #    -> we have to create a new list '[]' and insert the read_id/transcript_name in it '[read_id]'.
        # 2) all other appearances of a key:
        #    -> we just append new value to the list of transcripts/reads.

        if name_of_the_reference != '*':
            # if this is the first time we are adding reads to 'name_of_the_reference' key of the dictionary:
            # (.get returns either the value (list) under key 'name_of_the_reference' or 0 if there is no such key)
            if mapping_to_transcripts.get(name_of_the_reference, 0) == 0:
                mapping_to_transcripts[name_of_the_reference] = [read_id]
            # if it's not the first time:
            else:
                mapping_to_transcripts[name_of_the_reference].append(read_id)

            if current_split_dict_paf.get(read_id, 0) == 0:
                current_split_dict_paf[read_id] = [name_of_the_reference]
            else:
                current_split_dict_paf[read_id].append(name_of_the_reference)

        line = file_paf.readline()

    # saving dictionaries to pickle directory is the pafe in both cases so we will do it later on

else:

    # name of the file we will get our previous split's dictionary from
    pickle_file_name = ''

    # load previous split's dictionary
    if which_split == '2nd_split':
        pickle_file_name = '1st_split/read_mappings.p'

    elif which_split == '3rd_split':
        pickle_file_name = '2nd_split/read_mappings.p'

    elif which_split == '4th_split':
        pickle_file_name = '3rd_split/read_mappings.p'

    else:
        pickle_file_name = split_number[0] + 'th_split/read_mappings.p'

    with open(args.pickle_directory + '/' + pickle_file_name, 'rb') as f:
        previous_split_dict_paf = pickle.load(f)

    while line != '':
        parts = re.split('\s+', line)

        # parts[0] = id of the mapped read
        # parts[5] = name of the reference the read mapped to

        read_id = parts[0]
        name_of_the_reference = parts[5]

        if name_of_the_reference != '*':
            # After the first mapping, we are mapping unmapped reads to mapped ones. On account of that,
            # references in .paf file will be id's of previously mapped reads. To determine to which transcript bin
            # our read belongs ( or to witch transcript it 'maps' to ), we will use the previous split's dictionary.
            # From it, we can get the transcripts to which previously mapped read maps. Then we can say our read maps
            # to the pafe transcripts.

            transcripts_to_which_the_reference_mapped = previous_split_dict_paf[name_of_the_reference]

            for transcript in transcripts_to_which_the_reference_mapped:
                if mapping_to_transcripts.get(transcript, 0) == 0:
                    mapping_to_transcripts[transcript] = [read_id]
                else:
                    mapping_to_transcripts[transcript].append(read_id)

                if current_split_dict_paf.get(read_id, 0) == 0:
                    current_split_dict_paf[read_id] = [transcript]
                else:
                    current_split_dict_paf[read_id].append(transcript)

        line = file_paf.readline()

# saving dictionaries to pickle directory

# create current split's directory first

if not os.path.exists(args.pickle_directory + '/' + which_split):
    try:
        os.makedirs(args.pickle_directory + '/' + which_split, 0o700)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

# save files

with open(args.pickle_directory + '/' + which_split + '/transcript_bins.p', 'wb') as f:
    pickle.dump(mapping_to_transcripts, f)
with open(args.pickle_directory + '/' + which_split + '/read_mappings.p', 'wb') as f:
    pickle.dump(current_split_dict_paf, f)

# close PAF file

file_paf.close()

# read FASTQ/FASTA file and split mapped/unmapped reads

reads_to_split_file = open(args.reads_to_split, 'r', encoding='utf-8')
out_file_only_mapped = open(out_path_only_mapped, 'w', encoding='utf-8')
out_file_only_unmapped = open(out_path_only_unmapped, 'w', encoding='utf-8')


# if it is FASTQ file

if extension == '.fastq':

    # each entry in FASTQ file has 4 lines
    line1 = reads_to_split_file.readline()
    line2 = reads_to_split_file.readline()
    line3 = reads_to_split_file.readline()
    line4 = reads_to_split_file.readline()

    while line1 != '':
        parts = re.split('\s+', line1)

        read_id = parts[0][1:]

        # if read's id is in dictionary 'current_split_dict_paf' that means it mapped
        if current_split_dict_paf.get(read_id, 0) != 0:
            out_file_only_mapped.write(line1)
            out_file_only_mapped.write(line2)
            out_file_only_mapped.write(line3)
            out_file_only_mapped.write(line4)
        else:
            out_file_only_unmapped.write(line1)
            out_file_only_unmapped.write(line2)
            out_file_only_unmapped.write(line3)
            out_file_only_unmapped.write(line4)

        line1 = reads_to_split_file.readline()
        line2 = reads_to_split_file.readline()
        line3 = reads_to_split_file.readline()
        line4 = reads_to_split_file.readline()

# if it's FASTA file

elif extension == '.fasta':
    # each entry in FASTA file has 2 lines
    line1 = reads_to_split_file.readline()
    line2 = reads_to_split_file.readline()

    while line1 != '':
        parts = re.split('\s+', line1)

        read_id = parts[0][1:]

        # if read's id is in dictionary 'current_split_dict_paf' that means it mapped
        if current_split_dict_paf.get(read_id, 0) != 0:
            out_file_only_mapped.write(line1)
            out_file_only_mapped.write(line2)

        else:
            out_file_only_unmapped.write(line1)
            out_file_only_unmapped.write(line2)

        line1 = reads_to_split_file.readline()
        line2 = reads_to_split_file.readline()

reads_to_split_file.close()
out_file_only_mapped.close()
out_file_only_unmapped.close()
