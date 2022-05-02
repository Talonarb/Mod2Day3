# This section does the following:
# Imports necessary modules and packages
# The following two actions have defaults upon no user entry based on the location of the script.
# Takes user inputs for the directory of FASTA file to be assembled
# Takes user input for the directory of a query file 
# The following two actions have defaults upon no user entry, filenames default to those included in the repository
# Takes user input for the filename of a FASTA file to be assembled
# Takes user input for the filename of a query file to be used to find matches within the assembled contigs
Testing_multi_kmer = ["AAACTG", "CTGAAT", "AATCGG", "CGGAAAT"]
#This function takes the user inputs (or default values) and stores the information read from the files 
# Secondarily for some reason the QUERY file cannot be stripped of it's various newline and brackets when read 
import os
import itertools 
def Test_FASTAreader(ReadDirectory, QueryDirectory, ReadFilename, QueryFilename):
    global FASTAreads
    global FASTAquery
    ReadFile = open(os.path.join(ReadDirectory, ReadFilename))
    QueryFile = open(os.path.join(QueryDirectory, QueryFilename))
    FASTAreads = ReadFile.read()
    ReadFile.close
    FASTAquery = QueryFile.readlines()[1:]
    QueryFile.close
# This converts the line read from the query file (after skipping the identifier) to a string
    FASTAquery = str(FASTAquery)
# This step removes artifacts from the query file such that no brackets, newline characters, or quotations remain.
    FASTAquery = FASTAquery[1:-4]  
    return FASTAreads, FASTAquery

# Modified FASTA to dict function, via A.J. Uppal's code available at :
# https://stackoverflow.com/questions/29333077/reading-a-fasta-file-format-into-python-dictionary 
def Test_FASTA_to_dict(fastafile):
    global fasta_dict
    fasta_dict = {}
    with open(fastafile) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta_dict:
                    fasta_dict[active_sequence_name] = []
                continue
            sequence = line
            fasta_dict[active_sequence_name].append(sequence)
        #fastafile.close
    return(fasta_dict)

def fasta_dict_to_kv_lists(dictionary):
    global sequence_names_from_fasta
    global sequences_from_fasta
    sequence_names_from_fasta = list(dictionary.keys())
    sequences_from_fasta = list(itertools.chain.from_iterable(list(dictionary.values())))
    return(sequence_names_from_fasta, sequences_from_fasta)


def Variables_to_files(Name_of_variables_for_AllelesFasta, Name_of_Variable_for_AlleleALN):
    pass
# Variables_to_files(Name_of_variables_for_AllelesFasta, Name_of_Variable_for_AlleleALN)