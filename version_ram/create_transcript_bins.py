# To find out what this script does run the script with argument '-h' or '--help'.

from __future__ import print_function
import sys
import pickle
import os
import re
import shutil
import argparse
import errno

# function that allows printing to stderr


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# Functions that are used to check input formats and make output directories if needed:

def pickle_dir(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError('Pickle directory needs to exists and have all necessary files in it.')

    return x

def fa_or_fq_file(x):
    extensions = ('.fastq', '.fq', '.fasta', '.fa')

    if not x.lower().endswith(extensions):
        raise argparse.ArgumentTypeError("'reads_to_split' argument needs to be FASTA or FASTQ file.")

    return x


def output_dir(x):
    if not os.path.exists(x):
        try:
            os.makedirs(x, 0o700)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise

    return x

# function that checks if threshold is >=0


def threshold_type(x):
    try:
        x = int(x)
    except TypeError:
        raise argparse.ArgumentTypeError("Threshold needs to be an integer.")
    if x < 0:
        raise argparse.ArgumentTypeError("Minimum threshold is 0.")
    return x

# define arguments

# Argparse is a module that allows us to parse command line arguments easier:
# -> quick intro: https://docs.python.org/3/howto/argparse.html


parser = argparse.ArgumentParser(
    description=
    """
    This script merges dictionaries contained in all transcript_bins.p pickle files. Afterwards, it uses 
    the merged dictionary to split the reads we are mapping to transcriptome into transcript bins.
    
    Pickle files will contain dictionaries of following type:
        ### transcript_bins.p:
            transcript -> list of reads that mapped to it
        ### read_mappings.p: (will not be used)
            read -> list of transcripts it mapped to
    """

    , epilog=
    """
    Output:
        FASTQ/FASTA files named after transcripts which will contain all the reads that mapped to the 
        respective transcripts
    """
)
parser.add_argument(
    'pickle_directory', type=pickle_dir,

    help=
    """
    Path to the pickle files directory.
            
    NOTE:
       This argument is just the top directory of pickle files (The script requires directories 
       1st_split, 2nd_split, ... , <last_split>th_split to be present in that directory and will
       extract dictionaries from pickle files in each sub directory. These sub directories will
       be created by split_mapped_and_unmapped_reads.py script).
    """
)

parser.add_argument(
    'reads_to_split', type=fa_or_fq_file,

    help=
    """
    FASTQ/FASTA file containing all the reads. These reads will be split into bins according to the 
    transcript they mapped to. 
    
    NOTE: 
        Script will check if the file ends with '.fastq', '.fasta', '.fq', '.fa' or any of its variations.
        If not 'ArgumentTypeError' will be raised.
    """
)

parser.add_argument(
    'output_directory_for_transcript_bins', type=output_dir,

    help=
    """
    After running the script this directory will contain FASTQ/FASTA files with transcript names.
    Each file will harbour all the reads that mapped to the respective transcript. 
    If this directory doesn't exist the script will try to create it.
    """
)

parser.add_argument(
    '-t', '--threshold', type=threshold_type, default=0,

    help=
    """
    If the number of reads that mapped to a transcript is lower than this threshold, that transcript will 
    be ignored. If omitted or set to 0, then no threshold is used.
    """
)

parser.add_argument(
    '-d', '--delete_pickle_files', action='store_true',

    help=
    """
    If set, the script will, after splitting the reads into transcript bins, delete pickle files directory.
    """
)

### get arguments

args = parser.parse_args()


### merge dictionaries

# list in which we will store all dictionaries
dicts_transcript_bins = []

# merged dictionary
full_dict_transcript_bins = {}


### merge transcript_bins.p dictionaries

# get the dictionaries
for directory in os.listdir(args.pickle_directory):
    with open(args.pickle_directory + '/' + directory + '/transcript_bins.p', 'rb') as f:
        dicts_transcript_bins.append(pickle.load(f))

# we will merge all dictionaries into the first dictionary (pop removes the first dictionary and puts it in full_dict)
full_dict_transcript_bins = dicts_transcript_bins.pop()

# iterate through remaining dictionaries and merge them
for dictionary in dicts_transcript_bins:
    for (key, value) in dictionary.items():
        # if we have no key 'key' in full_dictionary we just put the list 'value' as the new key's value
        if full_dict_transcript_bins.get(key, 0) == 0:
            full_dict_transcript_bins[key] = value
        # else we just extend the value (which is a list) of key 'key' with list 'value'
        # extend will add all elements of one list to another list
        else:
            full_dict_transcript_bins[key].extend(value)


### remove transcripts from full_dict_transcript_bins with less then threshold reads mapped to them

if args.threshold > 0:
    for (key, value) in full_dict_transcript_bins.items():
        if len(value) < args.threshold:
            del full_dict_transcript_bins[key]


### split reads into transcript bins

# We will read the file 'reads_to_split' and make a dictionary of type:
#   read_id -> lines referring to the read in FASTQ or FASTA 'reads_to_split' file
#
# Then we will go through transcripts in 'full_dict_transcript_bins' dictionary. For each transcript we will create a
# new file in output directory named after the transcript. For each read that mapped to the transcript we will write
# the read's lines to the file.

dict_reads_to_split = {}

# first we check if reads_to_split is FASTQ or FASTA file
# then we split the reads

FAorFQ = ''

parts = re.split('/', args.reads_to_split)
parts2 = re.split('\.', parts[-1])
FAorFQ = parts2[1]

already_written = {}

# if it's .fasta or .fa
if FAorFQ == 'fasta' or FAorFQ == 'fa':
    with open(args.reads_to_split, 'r', encoding='utf-8') as input_reads:
        line1 = input_reads.readline()
        line2 = input_reads.readline()

        while line1 != '':
            # getting the read id
            parts = re.split('\s+', line1)
            read_id = parts[0][1:]

            dict_reads_to_split[read_id] = (line1, line2)

            line1 = input_reads.readline()
            line2 = input_reads.readline()

    # splitting the reads

    for (transcript, read_ids) in full_dict_transcript_bins.items():
            with open(args.output_directory_for_transcript_bins + '/' + transcript + '.fasta', 'w', encoding='utf-8') \
                    as output_file:
                for read_id in read_ids:
                    if read_id in already_written:
                        continue
                    
                    already_written[read_id] = True
                    lines = dict_reads_to_split[read_id]

                    output_file.write(lines[0])
                    output_file.write(lines[1])


# if it's .fastq or .fq
if FAorFQ == 'fastq' or FAorFQ == 'fq':
    with open(args.reads_to_split, 'r', encoding='utf-8') as input_reads:
        line1 = input_reads.readline()
        line2 = input_reads.readline()
        line3 = input_reads.readline()
        line4 = input_reads.readline()

        while line1 != '':
            # getting the read id
            parts = re.split('\s+', line1)
            read_id = parts[0][1:]

            dict_reads_to_split[read_id] = (line1, line2, line3, line4)

            line1 = input_reads.readline()
            line2 = input_reads.readline()
            line3 = input_reads.readline()
            line4 = input_reads.readline()

    # splitting the reads

    for (transcript, read_ids) in full_dict_transcript_bins.items():
        with open(args.output_directory_for_transcript_bins + '/' + transcript + '.fastq', 'w', encoding='utf-8') \
                as output_file:
            for read_id in read_ids:
                if read_id in already_written:
                        continue
                    
                already_written[read_id] = True

                lines = dict_reads_to_split[read_id]

                output_file.write(lines[0])
                output_file.write(lines[1])
                output_file.write(lines[2])
                output_file.write(lines[3])


### if delete_pickle_files is set, then delete pickle directory

if args.delete_pickle_files:
    try:
        shutil.rmtree(args.pickle_directory)
    except:
        eprint("ERROR: Couldn't delete pickle directory!")
        # exit(1)

### if delete_pickle_files not set, then save 'full_dict_transcript_bins' to pickle directory

else:
    with open(args.pickle_directory + '/transcript_bins_full.p', 'wb') as f:
        pickle.dump(full_dict_transcript_bins, f)
