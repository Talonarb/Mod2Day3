# Create main function, call the various functions in order, such that the variables are passable to one another without being global
# aside from the main reads. 

# This function assumes that the incoming FASTA files are small enough to fit in memory and are portable. This
# is not designed for files larger than ~500 MB

# Add to readme: if you want to visualize the graph use the python jupyter notebook included, and have the
# extensions gvmagic and graphviz. Secondarily, for larger sequences it is untenable to visualize the graph/
# tree, as it would not be human readable. Exact matching used, is assumed that the sequence of interest must
# be found in its entirety. The matching portion assumes all sequences have only the alphabet of ACTG. 

# Importing built-ins

import os.path     
import pandas as pd 
import sys
import re 
import numpy as np
import bisect
import string
from itertools import permutations
sys.path.insert(0, (os.getcwd()))

# Importing scripts
from src import DependenciesAndFileImport_TEST as DandF, DuhBrewin as Dbrew, Heirholzer_test as Hhzer, BoyerMoore as BoyM 

# This function was having trouble being imported properly, so it defined within the main script. 

# Edge_dict is a generated list of edges from the edge list of the kmerizaiton function
# I want to call the dictionary generation based on the input to the pathfinding algorithm, and return a dictionary for use within
# the pathfinding function.
# Convert list of edges from DeBruijn function into a dictionary of source node as key and value pairs as directed targets of nodes
def pairs_to_dict(pair_list, dictionary_to_be_made): 
    for a, b in pair_list: # loops through list of pairs 
        dictionary_to_be_made.setdefault(a, []).append(b) # makes a key of the "a" position in pairs with value "b" 
    return dictionary_to_be_made # and now we have our dictionary



# This is where global variables and the working directory are set. Please ensure all files are in the current working directory
# for importing ths scripts, although the directory of the READS.fasta and QUERY.fasta files can be specified as elsewhere. 
os.chdir(os.path.dirname(sys.argv[0]))
ReadDirectory = (input("Enter the directory of the FASTA reads file: ")) or (os.path.join(os.getcwd()))
QueryDirectory = (input("Enter the directory of the QUERY query file: ")) or (os.path.join(os.getcwd()))
ReadFilename = (input("Enter the filename of the FASTA reads file: ")) or ("READS.fasta")
QueryFilename = (input("Enter the filename of the  reads file: ")) or ("QUERY.fasta")
ReadFile = open(os.path.join(ReadDirectory, ReadFilename))
QueryFile = open(os.path.join(QueryDirectory, QueryFilename))
FASTAreads = ReadFile.read()
ReadFile.close
FASTAquery = QueryFile.readlines()[1:]
QueryFile.close 

#Kmer size selection
kmer_size = (input("Enter the desired Kmer size, as an integer. The default is 4, which is what is used for the test file : ")) or (4)
# And convert the user input from string format into int, or ensure that the defaul is read as int as opposed to string. 
kmer_size = int(kmer_size)

# This is for graph traversal, as starting nodes impact the algorithm quite a lot, thus users can change which starting node should be used. In future work it may be advisable to 
# find a list of potential starting nodes using incoming/outgoing edges, such that nodes with no outgoing edges are not selectable, and nodes with only outgoing edges are preferred
# and a list is provided to the user to choose from. 
starting_node = (input("Enter the desired starting node, as an integer. The default is 1, which is what is used for the test file : ")) or (1)
starting_node = int(starting_node)
# uncomment this section of the code if you desire to run the pipeline on the testfile included in this repository. 

#TestSearchFilename = "TESTINGSEARCH.fasta"
#TestSearch = open(os.path.join(ReadDirectory, TestSearchFilename))
#TestSearchReads = TestSearch.readlines()


# This converts the line read from the query file (after skipping the identifier) to a string
FASTAquery = str(FASTAquery)
print('This is the fastquery before artifact removal, ', FASTAquery)
# This step removes artifacts from the query file such that no brackets, newline characters, or quotations remain. 
# Alter this section if there are no such characters remaining in your query fasta file. 
print(QueryFilename,'', ReadFilename)
print(QueryDirectory,'', ReadDirectory)
def main():

    import os
    global circuit 
    FASTAreads, FASTAquery = DandF.Test_FASTAreader(ReadDirectory, QueryDirectory, ReadFilename, QueryFilename)
    FASTAquery = FASTAquery[2:-4]
    print("this is the fastaquery ", FASTAquery)
    FASTAquery_reverse = FASTAquery[::-1] 
    fasta_dict = DandF.Test_FASTA_to_dict((os.path.join(ReadDirectory, ReadFilename)))
    sequence_names_from_fasta, sequences_from_fasta = DandF.fasta_dict_to_kv_lists(fasta_dict)
   
    # uncomment this section of the code if you desire to run the pipeline on the testfile included in this repository. 
    # test_edge_dict = {}
    # test_sequence_dict = DandF.Test_FASTA_to_dict((os.path.join(ReadDirectory, TestSearchFilename)))
    #test_nodes, test_edges = Dbrew.DuhBrewin(test_seq, kmer_size)
    #Dbrew.graph_representation(test_seq, kmer_size)
    #pairs_to_dict(test_edges, test_edge_dict)
    #testcontig = Hhzer.visualize_Eulerian_tour(test_edges, test_edge_dict)
    #print(testcontig)

    # This section creates the Kmers, DeBruijn graph, and creates a contig. Several sizes of K as well as 
    # starting nodes may be necessary for achieving the largest possible contig. 
    FASTA_edge_dict = {}
    FASTA_nodes, FASTA_edges = Dbrew.DuhBrewin(sequences_from_fasta, kmer_size)
    pairs_to_dict(FASTA_edges, FASTA_edge_dict)
    FASTA_contig = Hhzer.visualize_Eulerian_tour(FASTA_edges, FASTA_edge_dict, starting_node)
    print("This is the assembled contig, ",  FASTA_contig)
    print("this is the fastaquery, ", FASTAquery)
    p_bm = BoyM.BoyerMoore(FASTAquery)
    Match_index = (BoyM.boyer_moore(FASTAquery, p_bm, FASTA_contig))
    if Match_index == []:
        print("The query is not found in the contigs ")
    else: print("The index of the start of the query match to the contig is ", Match_index)

    p_bm = BoyM.BoyerMoore(FASTAquery_reverse)
    Match_index = (BoyM.boyer_moore(FASTAquery_reverse, p_bm, FASTA_contig))
    if Match_index == []:
        print("The reveresed query is not found in the contigs ")
    else: print("The index of the start of the reversed query match to the contig is ", Match_index)

    d = {'SSEQID': [], 'QSEQID': [], 'SSTART': [], 'SEND': [], 'QSTART': [], 'QEND': []}
    outdataframe = pd.DataFrame(data = d)
    print(outdataframe)
    
    # This would be another loop for each contig if there are more than one contigs built, nested for loop
    # when appending contig, you append the contig_id for each contig being looped through, or contig + [i/j], whatever I loop through
    # set contig count to 0, contig count +=1 
    for sequence_id_from_dict, subject_sequence in fasta_dict.items():
        # if no match, skip 
        converted_seq = '' .join(subject_sequence)
        p = converted_seq
        p_bm = BoyM.BoyerMoore(p)
        Match_index = (BoyM.boyer_moore(p, p_bm, FASTA_contig))
        if Match_index != []:
            # record sequence id in SSEQID
            # record contig ID in QSEQID contig_1 in this case
            # record SSTART index for where the sequence begins to align with respect to the sequence
            # record SEND for where the sequence alignment ends with respect to the sequence
            # QSTART as above, but with resepect to the contig being aligned to
            # QEND as above, but with resepect to the contig being aligned to
            match_index_to_int = Match_index[0]
            data_to_append = [sequence_id_from_dict, 'contig_found', '0', len(converted_seq), match_index_to_int, ((len(converted_seq) + match_index_to_int))]
            outdataframe = outdataframe.append(pd.DataFrame([data_to_append], columns = ['SSEQID', 'QSEQID', 'SSTART', 'SEND', 'QSTART', 'QEND']), ignore_index = True)
            print(data_to_append)
            print(outdataframe)
            outdataframe.to_csv('ALLELES.aln', sep="\t")

    



if __name__ == "__main__" :
    main()