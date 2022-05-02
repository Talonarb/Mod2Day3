# The implementation of alignment is now focused on the alignment itself, using K-mer compositions of the query sequence and the sequence
# to be aligned to. K-mer composition is then indexed after finding all unique K-mer compositions for both targets, and building a 
# reference hashtable for both for alignment purposes. As the goal is only the alignment of the longest contig that contains the query sequence
# this can be sped along as there is a fair chance some reads will not be involved in the final assembled contig. 

# This also doubles as the alignment indexer for each read, which can be stored according to the read identity in the FASTA file which
# is stored in a dictionary in the first steps of this program. The hashtable can then be queried using a binary search, making lookup of
# where Kmers are located rapid, assuming of course after the initial hashtable creation the data is sorted by Kmer alphabetically, which
# would then immediately split search space in half for every new search it conducts until the hashtable is exhausted or the kmer location
# is found. The hashtable is built by first querying the unique K-mers of a given size for a sequence (instead of computing all possible
# compositions of Kmers only those present in the sequence are used), following the sliding window approach of finding all unique Kmers, 
# those kmers are used to create keys for a dictionary, and while that happens the dictionary updates the value pairs for each key by
# referencing the location that the Kmer occurs. 

# This assumes that the kmers are more or less unique to a given sequence. Could be improved. This would not be good for a highly 
# repetitive sequence. Set assumptions in: readme and report. 
import bisect 
class Index(object):
    # Input t is text, k is the value of k to be used for kmer composition creation
    def __init__(self, t, k):
        self.k = k
        self.index = []
# Provides every index for the sequence passed in, without going beyond the length of the sequence being indexed
        for i in range(len(t) - k +1):
            self.index.append((t[i:i+k], i))
        self.index.sort()
        return self.index 
# Without sorting binary search cannot be reasonably used for lookup, even if the sorting process may take a while it makes the exploration
# of the search space much faster as the search space is no longer O(n), but O(log base 2 n). 
# See: https://www.rudikershaw.com/articles/whichsearch 

# This function takes the index list created above as input as well as a kmer designated as "p" for pattern, and can then query the 
# sorted list to find the matches in a given string. 
    def query(self, p):
        kmer = p[:self.k] # This takes the pattern (p), and find the first k-bases of p to find where this pattern is in the index 
        i = bisect.bisect_left(self.index, (kmer, -1)) # This just finds the position this patterns occurs, assuming the positions is >-1
        kmer_matches = []
        while i < len(self.index):
# This break is for stopping if the location in the index is not the kmer, we can stop, as we know the index is sorted, so the matches are
# guaranteed to be exhausted after we pass the final occurence of our Kmer 
            if self.index[i][0] != kmer:
                return kmer_matches
# break is bad in programs, if breaks the program, it just stops everything. Use return in lieu of break. Don't do this in prelims. 

# This is the instance of where there is a match of the first k characters of our pattern p in to the text we are querying, and then 
# the position at which this match occurs is appended to the index. 
            kmer_matches.append(self.index[i][1])
            i += 1
        return kmer_matches

# This takes the first k units of a query pattern (p), matches it against a text to be queried (t), and then finds all the offsets 
# of those matches, records them into a list, and more importantly also checks that not only does the first k characters match the 
# text being queried, but also ensures that the remaning bases in the text being queried with pattern p match for the length of p. 
def query_Index(p, t, index):
    k = index.k
    offsets = []
    for i in index.query(p):
        if p[k:] == t[ i + k :i + len(p) ]:
            offsets.append(i)
    return offsets