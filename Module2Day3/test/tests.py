# Unit tests for all the various functions. Inspired by help from Katerina Cortes, Erik Serrano, Sofia Colorado,
# Nicholas Garcia, and with inspiration via Dr. Ryan Layer (tests are good, check out his courses). Importing modules.

import unittest
import os.path     
import sys
import re 
import bisect
import string   
import itertools 

# Sets path to current directory, then the parent directory of current directory. 
current_direc = (os.path.dirname(sys.argv[0]))
print(current_direc)
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
print(os.getcwd())
sys.path.insert(0, (os.getcwd()))

# Importing scripts 
from src import DependenciesAndFileImport_TEST as DandF, DuhBrewin as Dbrew, Heirholzer_test as Hhzer, BoyerMoore as BoyM 
# import DependenciesAndFileImport_TEST as DandF 
# import DuhBrewin as Dbrew
# import Heirholzer_test as Hhzer
# import BoyerMoore as BoyM


# This is where global variables and the working directory are set. Please ensure all files are in the current working directory
# for importing ths scripts, although the directory of the READS.fasta and QUERY.fasta files can be specified as elsewhere. 
ReadDirectory = (input("Enter the directory of the FASTA reads file: ")) or (os.path.join(os.getcwd()))
TestSearchFilename = "TESTINGSEARCH.fasta"
TestSearch = open(os.path.join(ReadDirectory, TestSearchFilename))
TestSearchReads = TestSearch.readlines()
TestSearch.close 

def pairs_to_dict(pair_list, dictionary_to_be_made): 
    for a, b in pair_list: # loops through list of pairs 
        dictionary_to_be_made.setdefault(a, []).append(b) # makes a key of the "a" position in pairs with value "b" 
    return dictionary_to_be_made # and now we have our dictionary

# First test is for file reading and dictionary creation. Simple task, but important to check. 
# Tests both the file reading and dictionary creation functions. 
class FastaReaderTest(unittest.TestCase):
    
    def test_reader(self):
        reads = DandF.Test_FASTA_to_dict((os.path.join(ReadDirectory, TestSearchFilename)))
        self.assertEqual(reads, {'TEST2:FASTAFILE:TEST': ['TATTGTT'], 'TEST3:FASTAFILE:TEST': ['GCATGG'], 'TEST4:FASTAFILE:TEST': ['GCATGCA'], 'TEST5:FASTAFILE:TEST': ['TGGCTC'], 'TEST7:FASTAFILE:TEST': ['ATTGTTTTTAGA']})

# Second test is to make sure the dictionary is being properly converted to lists with matching indexes for
# each key/value pair. This is used in later functions. Uses index matching, so important they line up. 
class DictToList(unittest.TestCase):
    
    def test_DictToList(self):
        reads = DandF.Test_FASTA_to_dict((os.path.join(ReadDirectory, TestSearchFilename)))
        sequence_names_from_test, sequences_from_test = DandF.fasta_dict_to_kv_lists(reads)
        self.assertEqual(sequence_names_from_test, ['TEST2:FASTAFILE:TEST', 'TEST3:FASTAFILE:TEST', 'TEST4:FASTAFILE:TEST', 'TEST5:FASTAFILE:TEST', 'TEST7:FASTAFILE:TEST'])
        self.assertEqual(sequences_from_test, ['TATTGTT', 'GCATGG', 'GCATGCA', 'TGGCTC', 'ATTGTTTTTAGA'])

# Next steps in the pipeline are kmer composition creation, DeBruijn graph visualization, and traversal. 
class DuhBrewin(unittest.TestCase):
    
    def test_kmerization(self):
        sequences_from_test = ['TATTGTT', 'GCATGG', 'GCATGCA', 'TGGCTC', 'ATTGTTTTTAGA']
        nodes, edges = Dbrew.DuhBrewin(sequences_from_test, 4)
        self.assertEqual(nodes, {'GCA', 'TAG', 'ATG', 'TGT', 'CTC', 'TGC', 'AGA', 'GTT', 'GCT', 'TGG'})
        self.assertEqual(edges, [('TAT', 'ATT'), ('ATT', 'TTG'), ('TTG', 'TGT'), ('TGT', 'GTT'), ('GCA', 'CAT'), ('CAT', 'ATG'), ('ATG', 'TGG'), ('GCA', 'CAT'), ('CAT', 'ATG'), ('ATG', 'TGC'), ('TGC', 'GCA'), ('TGG', 'GGC'), ('GGC', 'GCT'), ('GCT', 'CTC'), ('ATT', 'TTG'), ('TTG', 'TGT'), ('TGT', 'GTT'), ('GTT', 'TTT'), ('TTT', 'TTT'), ('TTT', 'TTT'), ('TTT', 'TTA'), ('TTA', 'TAG'), ('TAG', 'AGA')])

class HeirHolzer(unittest.TestCase):
    
    # This test is designed to make sure the path chosen for a given set of inputs always produces the same contig so long as the 
    # starting node and edge list/dictionary are the same. The original sequence was  "TATTGT", and k-1mers are desgined to have
    # no duplicates or other paths, such that it will always assembly the same contig. This is not always true when k-1mers have
    # several repetitive elements, or in the case of many repetitions of a k-1mer. Contig assembly varies in those conditions, and
    # the output contig is not necessarily always the same. 
    def test_heirholzer(self):
        test_unit_edge_dict = {}
        testing_stuff_edges = [('TAT', 'ATT'), ('ATT', 'TTG'), ('TTG', 'TGT')]
        pairs_to_dict(testing_stuff_edges, test_unit_edge_dict)
        contig_test_check = Hhzer.visualize_Eulerian_tour(testing_stuff_edges, test_unit_edge_dict, 1)
        self.assertEqual(contig_test_check, 'TATTGT')

class BoymerMoore(unittest.TestCase):

    # This test is a simple matching test to find out whether the Boyer-Moore implementation for string matching is working correctly. 
    # This function relies on a series of calculated arrays to be used before hand, and is predicated on the idea that only the 
    # characters "A", "C", "T", and "G" are potential in the strings being compared. The index the pattern matches the string to be
    # compared to should be "8" and "17", as 8 is the index of the beggining of the first occurence of "TGGCTA" in the string, and 17 
    # is the beginning of the next occurence in the string. Note that this occurences are stored as a list, such that several 
    # occurences can be found and recorded. 
        def test_Boyer_Moore(self):
            contig_found = "TCCTGGCTTGGCTACTTTGGCTA"
            pattern = "TGGCTA"
            p_bm = BoyM.BoyerMoore(pattern)
            occurence_of_match = (BoyM.boyer_moore(pattern, p_bm, contig_found))
            self.assertEqual(occurence_of_match, [8,17])
        





if __name__ == '__main__':
    unittest.main()
