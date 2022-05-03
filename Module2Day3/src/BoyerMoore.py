# Boyer moore matching algorithm, inspired by ADS1: Practical: implementing Boyer-Moore on youtube, accesible at the following link:
# https://youtu.be/CT1lQN73UMs 
# Note that this is a playlist, and has several components 

# Pre-processes the reads to find the positional arguments of where the text mismatches, and either a mismatch of an alignment becomes
# a match, or the character index of the mismatch is passed. This is the "bad character" rule, and when combined with the "good suffix"
# rule, makes alignment faster. The "good suffix" rule uses a subset of the string being aligned that is repeated within the full string, and
# aligns using the next matching subset within the string. This will continue as the suffix continues to improve alignment stepwise, 
# until either the alignment fully matches or all possible alignments are exhausted, and no complete alignment is found. When combined, 
# the calculated values of shifts for either rule is used to choose which to proceed with, and whichever number of alignments skipped
# is larger is used. Hence, many improper alignments will be skipped because they cannot be correct. This allows characters that are not
# part of the text to be aligned to be skipped, as well as knowably incorrect alignments to be skipped. Secondarily, the number of skips
# can be precalculated using the pattern that is being aligned. Thus, the choice of whether to use either rule can be predetermined, and 
# a lookup table can be created to further speed up the comparison process, by using a table of length n, which corresponds to the length
# of the query string (not the string of characters it is being aligned to) by 4, the 4 being the possible nucleotides it can mismatch
# for when it is being aligned. This could be implemented for other uses, e.g., protein alignment, and instead have an n x 20 table. 

## A theorem called the "Gusfield theorem", the principals of which explain how to do preprocessing of alignment
## while using the Boyer-Moore algorithm, is explained in detail at the following URL: 
## http://www.cbcb.umd.edu/confcour/Spring2010/CMSC858W-materials/restrict/Gusfield-0-1-2-3.3.pdf
## I will try to use terms described above throughout to show what each section of the code is doing, as each part is inspired by the text
## linked above as well as the youtube linked above. 

# This is the really neato pre-processing part, where the initial values of a string to be used for alignment can be calculated for the
# "good suffix" or "bad character" rules, using a bunch of stuff from the Gusfield theorem, which basically involved calculating a series
# of arrays and leveraging them to make the lookup table for when to use which rule. 

# This approach may be less appropriate than kmer indexing and tracking, but currently that is something I do not yet know how to implement,
# depending on how the meeting with my tutor goes I may be able to take that approach instead, which seems more efficient. 

# This creates the z array from a given string (s), which is used in the creation of a gang of other arrays which eventually result in a 
# table of when to implement uses of each rule based on mismatch information. 


from ast import Assert


def z_array(s):

    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    for i in range( 1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else: 
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range (2, len(s)):
        assert z[k] == 0 
        if k > r:
            for i in range (k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            nbeta = r - k + 1
            zkp = z[k -1 ]
            if nbeta > zkp:
                z[k] = zkp
            else: 
                nmatch = 0
                for i in range (r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z

# Compiles the N array (Gusfield Theorem) from the Z array

def n_array(s):
    return z_array(s[::-1][::-1])

def big_l_prime_array(p, n):

    # Compile L prime array (gusfield theorem) using p and N array. 
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp

def big_l_array(p, lp):

# Compiles the L array described in the Gusfield theorem, using p and L prime array.
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):

    # Compile lp' array (gusfield theorem 2.2.4) using N array
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i + 1:
            small_lp[len(n)-i-1] = i + 1 ## this is the prefix/suffix matching aspect, which is part of the "good suffix" rule
    for i in range(len(n)-2, -1, -1): 
        if small_lp[i] == 0:
            small_lp[i] == small_lp[i + 1]
    return small_lp

def good_suffix_table(p):

    # Gives us the tables that we need to know when to apply the good suffix rule #
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)

def good_suffix_mismatch(i, big_l_prime, small_l_prime):
# Using L array, L' array, l array, and l' array, gives the amount to shift when using the good suffix rule
    length = len(big_l_prime)
    assert i < length
    if i == length -1:
        return 0
    i += 1 
    # i in this scenario is the patterns earliest match in comparison to the sequence being aligned to
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]

def good_suffix_match(small_l_prime):
# In the instance of a complete match of the query sequence to the sequence to be aligned to, output the shift calculated based on the 
# "good suffix" rule. 
    return len(small_l_prime) - small_l_prime[1]

def dense_bad_char_tab(p, amap):
# Using the query sequence and a list with the possible characters (sorted), create and return a table of the "bad characters".
# The table is indexed by offset, then by character
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab 

class BoyerMoore(object):
# Encapsulates pattern and associated Boyer-Moore preprocessing

    def __init__(self, p, alphabet = 'ACGT'):
        self.p = p
        self.alphabet = alphabet

#Makes a map from the alphabet to integers

        self.amap = {}
        for i in range(len(self.alphabet)):
            self.amap[self.alphabet[i]] = i
        
 # Bad character table creation
        self.bad_char = dense_bad_char_tab(p, self.amap)

# Good suffix table creation
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i , c):
# Outputs the the number of skips using the bad character rule at a given offset, i # 
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci]-1)
        return i - (self.bad_char[i][ci]-1)
    
    def good_suffix_rule(self, i):
# Given a mismatch at offest i, return amount to shift as determed by the good suffix rule
        length = len(self.big_l)
        assert i < length
        if i == length -1:
            return 0
        i += 1 # Again, earliest match of a given query sequence to a sequence to be compared to in this case
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length -self.small_l_prime[i]
    
    def match_skip(self):
# Outputs the amount to shift when there is an exact match
        return len(self.small_l_prime) - self.small_l_prime[1]


# This function does pattern matching based on bad character alignments or good suffix alignments, with arguments of the pattern, the
# pattern's precalculated boyermoore skip rule sizes, and t, being text to compare the pattern to. 
def boyer_moore(pattern, p_bm ,t):
    i = 0
    occurrences = []
    while i < len(t) - len(pattern) +1 :
        shift = 1
        mismatched = False
        for j in range(len(pattern)-1, -1, -1):
            if not pattern[j] == t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift,skip_gs)
        i += shift 
    return(occurrences)
