# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 20:35:13 2015
This project includes findPatternV2 class. The findPatternV2
class finds the occurance and position of a given pattern in a given genomic 
sequence.
@author: zhihuixie
"""
import bisect

class findPatternV2 ():
    """
    This class finds the occurance and position of a given pattern in a given 
    genomic sequence in a file.
    """
    def __init__(self, pattern, filename = False):
        # initiate parameters
        self.pattern = pattern
        self.filename = filename
        
    def readGenome (self):
        """
        read genomic DNA sequence to a string
        """
        genome = ""
        with open (self.filename, "r") as f:
            for line in f:
                # skip header
                if not line[0] == ">":
                    genome += line.rstrip()
            f.close()
        return genome
        
    def naiveMatch (self, numberOfMismatch, text = False):
        """
        this is naive match to find the index of matched patterns in a genome
        and calculate number of total character comparisons and alignments
        """
        if text != False: #for test cases
            genome = text
        else:
            genome = self.readGenome()
        pattern = self.pattern
        occurances = []
        alignments = 0
        comparisons = 0
        for i in range(len(genome) - len(pattern) + 1):
            match = True
            counter = 0
            for j in range(len(pattern)):
                comparisons += 1
                if pattern[j] != genome[i+j]:
                    counter += 1
                if counter > numberOfMismatch:
                    match = False
                    break
            if match:
                occurances.append(i)
            alignments += 1
        return occurances, alignments, comparisons
        
    def boyerMoore (self, numberOfMismatch, bm, text = False):
        """
        this is naive match to find the index of matched patterns in a genome
        and calculate number of total character comparisons and alignments
        """
        i = 0
        if text != False: #for test cases
            genome = text
        else:
            genome = self.readGenome()
        pattern = self.pattern
        occurances = []
        alignments = 0
        comparisons = 0
        while i < len(genome) - len(pattern) + 1:
            shift = 1
            match = True
            for j in range(len(pattern) - 1, -1, -1):
                comparisons += 1
                if pattern[j] != genome[i+j]:
                    badCharacterSkip = bm.bad_character_rule(j, genome[i+j])
                    goodSuffixSkip = bm.good_suffix_rule(j)
                    shift = max(shift, badCharacterSkip, goodSuffixSkip)
                    match = False
                    break
            if match:
                occurances.append(i)
                goodSuffixSkip = bm.match_skip()
                shift = max(shift, goodSuffixSkip)
            i += shift
            alignments += 1
        return occurances, alignments, comparisons
        
    def matchedIndex(self, index, k_mer, pattern, isSubseqIndex):
        """
        find number of hits, occurances and time of occurance for a given pattern
        using string index in genome
        """
        genome = self.readGenome()
        occurances_match = []
        hit_index = []
        occurance_genome = []
        counter = 0
        if not isSubseqIndex:
            length = len(pattern)-k_mer + 1
        else:
            length = isSubseqIndex
        for i in range(length): # loop over to generate kmers
            if not isSubseqIndex:                
                pattern_q = pattern[i:i+k_mer]
            else:
                pattern_q = pattern[i:]
            hits = index.query(pattern_q) # query each kmer
            for hit in hits:
                counter += 1 #count total number of hits
                text = genome[hit-i : hit+len(pattern) -i]
                if hit-i not in hit_index: #avoid duplicated counts
                    hit_index.append(hit-i)
                    occurance, _, _ = self.naiveMatch(2,text)
                    occurances_match.extend(occurance)
                if len(occurance) != 0 and hit-i not in occurance_genome:
                    occurance_genome.append(hit-i) 
        return occurance_genome, len(occurances_match), counter
        
class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
    def genome_index(self):
         return self.index       

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

if __name__ == "__main__":
    from bm_preproc import BoyerMoore
    #Questions 1-3
    filename = ("../data/chr1.GRCh38.excerpt.fasta")
    #Q1: How many alignments does the naive exact matching algorithm try when 
    #matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG 
    #(derived from human Alu sequences) to the excerpt of human chromosome 1? 
    #(Don't consider reverse complements.)
    pattern = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
    patterns = findPatternV2(pattern, filename)
    print "Q1: The alignments for naive match algorithm is %d\n"%patterns.naiveMatch(0)[1]
    patterns = findPatternV2("GGCGCGGTGGCTCACGCCTGTAAT", filename)    
    
    #Q2: How many character comparisons does the naive exact matching algorithm 
    #try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG 
    #(derived from human Alu sequences) to the excerpt of human chromosome 1? 
    #(Don't consider reverse complements.)
    print "Q2: The characters comparisons for naive match algorithm is %d\n"%patterns.naiveMatch(0)[2]
    
    #How many alignments does Boyer-Moore try when matching the string 
    #GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG 
    #(derived from human Alu sequences) to the excerpt of human chromosome 1? 
    #(Don't consider reverse complements.)
    print "Q3: The alignments for Boyer-Moore algorithm is %d\n"%patterns.boyerMoore(0, \
          BoyerMoore(pattern, "ACGT"))[1]
    
    #Q4: How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, 
    #which is derived from a human Alu sequence, occur with up to 2 
    #substitutions in the excerpt of human chromosome 1? 
    #(Don't consider reverse complements here.)
    k_mer = 8      
    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    genome = patterns.readGenome()
    index = Index(genome, k_mer)
    occurances, numberOfOccurs, numberOfhits = patterns.matchedIndex(index, \
                                               k_mer, pattern, isSubseqIndex = False)
    print "Q4: Within 2 mismatchs, the string occurs %d times\n" %numberOfOccurs
    
    #Q5:Using the instructions given in Question 4, how many total index hits 
    #are there when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with 
    #up to 2 substitutions in the excerpt of human chromosome 1?
    print "Q5: Within 2 mismatchs, the total index hits are %d \n" %numberOfhits
    
    #Q6: When using this function, how many total index hits are there when 
    #searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the 
    #excerpt of human chromosome 1? (Again, don't consider reverse complements.)
    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    k_mer = 8
    vial = 3
    index = SubseqIndex(genome, k_mer, vial)
    occurances, numberOfOccurs, numberOfhits = patterns.matchedIndex(index, \
                                               k_mer, pattern, isSubseqIndex = vial)

    print "Q6: Within 2 mismatchs, the hits are %d\n" %numberOfhits
    