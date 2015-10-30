# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 19:55:08 2015
This project includes findPatternV3 class. The findPatternV3
class finds the edit distance of a given pattern in a given genomic 
sequence and constructs overlap graphs.
@author: zhihuixie
"""

from itertools import permutations

class findPatternV3 ():
    """
    This class finds the edit distance of a given pattern in a given genomic 
    sequence and constructs overlap graphs.
    """
    def __init__(self, filename, pattern = False):
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
        
    def readFastq(self):
        """
        read dna sequence and quality base from a fastq sequencing file to lists
        """
        with open (self.filename, "r") as f:
             sequences = []
             qualities = []
             while True:
                 f.readline() # skip name line
                 seq = f.readline().rstrip() # read sequence line
                 f.readline() # skip strand line
                 qual = f.readline().rstrip() # read quality line
                 if len(seq) == 0: #finish read
                     break
                 # add seqence and quality information to list
                 sequences.append(seq)
                 qualities.append(qual)
             f.close()
        return sequences, qualities
        
    def editDistance(self):
        """
        Implement dynamic algorithm to calculate edit distance between a given
        pattern and a given genome
        """
        pattern = self.pattern
        genome = self.readGenome()
        pattern_length = len(pattern) + 1
        genome_length = len(genome) + 1
        #generate matrix
        matrix = [[0]*genome_length for i in range(pattern_length)]
        # initiate the first column
        for i in range(pattern_length):
            matrix[i][0] = i
        for i in range(1, pattern_length):
            for j in range(1, genome_length):
                dist_hor = matrix[i][j-1] + 1
                dist_vel = matrix[i-1][j] + 1
                dist_diag = matrix[i-1][j-1] + 1 if pattern[i-1] != genome[j-1]\
                            else matrix[i-1][j-1]
                matrix[i][j] = min(dist_hor, dist_vel, dist_diag)
        return min(matrix[-1])                
            
    def phraseReads(self, k_mer):
        """
        construct the prefix and suffix of a read to a dictornary with read as
        key and pre-,suffix as values
        """
        reads, _ = self.readFastq()
        reads_dict = {}
        for read in reads:
            for i in range(len(read) - k_mer + 1):
                substring = read[i:i+k_mer]
                if substring not in reads_dict:
                    reads_dict[substring] = set([read])
                else:
                    reads_dict[substring].add(read)
        return reads_dict
    def overlap(self, read1, read2, k_mer):
        """
        find overlaped leftmost offset
        """
        start = 0
        while True:
            start = read1.find(read2[:k_mer], start)
            if start == -1:
                return 0 # without overlap
            if read2.startswith(read1[start:]):
                return len(read1) - start
            start += 1
            
    def overlapGraph(self, k_mer):
        """
        construct graph with key as a read (node) and values as all other 
        reads overlapped with the previous read (node)
        """
        reads_dict = self.phraseReads(k_mer)
        reads, _= self.readFastq()
        graph = {}
        for read1 in reads:
            k_mer_string = read1[len(read1) - k_mer:]
            if k_mer_string in reads_dict:
                edges = set([])
                reads_set = reads_dict[k_mer_string]
                for read2 in reads_set:
                    if read1 != read2: #skip self comparison
                        offset = self.overlap(read1, read2, k_mer)           
                        if offset > 0: # skip non-overlapped pairs
                            edges.add(read2) #add overlapped reads to be values
                            graph[read1] = edges           
        return graph
                        
    def naive_overlap_map(self, k_mer):
        """
        construct graph with key as a pair of reads with overlap and values as
        the leftmost offset of the overlap
        """
        graph = {}
        reads, _ = self.readFastq()
        for read1, read2 in permutations(reads, 2):
            #skip non-overlapped reads
            if read1[len(read1) - k_mer:] in read2:
                offset = self.overlap(read1, read2, k_mer)
                    # check if reads[i] overlapped with reads[j]
                if offset != 0:
                    graph[(read1,read2)] = offset  
        return graph

if __name__ == "__main__":
    
    #Q1: What is the edit distance of the best match between pattern 
    #GCTGATCGATCGTACG and the excerpt of human chromosome 1? 
    #(Don't consider reverse complements.)
    pattern = "GCTGATCGATCGTACG"
    filename = "../data/chr1.GRCh38.excerpt.fasta"
    patterns = findPatternV3 (filename, pattern)
    edit_dist = patterns.editDistance()
    print "Q1: the edit distance of the best match between pattern and the genome is %d\n"\
           %edit_dist
           
    #Q2: What is the edit distance of the best match between pattern 
    #GATTTACCAGATTGAG and the excerpt of human chromosome 1? 
    #(Don't consider reverse complements.)
    pattern = "GATTTACCAGATTGAG"
    filename = "../data/chr1.GRCh38.excerpt.fasta"
    patterns = findPatternV3 (filename, pattern)
    edit_dist = patterns.editDistance()
    print "Q2: the edit distance of the best match between pattern and the genome is %d\n"\
           %edit_dist
    
    #Q3: Picture the overlap graph corresponding to the overlaps just calculated. 
    # How many edges are in the graph? In other words, how many distinct pairs 
    # of reads overlap?
    #Q4: Picture the overlap graph corresponding to the overlaps computed for 
    #the previous question. How many nodes in this graph have at least one 
    #outgoing edge? (In other words, how many reads have a suffix involved in 
    #an overlap?)
    import time
    t1 = time.time()
    filename = "../data/ERR266411_1.for_asm.fastq"
    patterns = findPatternV3 (filename)
    k_mer = 30
    graph = patterns.naive_overlap_map(k_mer)
    t2 = time.time()
    print "Running time for naive overlap mapping: %d sec\n"%(t2 - t1)
    
    reads = patterns.phraseReads(k_mer)
    t3 = time.time()
    print "Running time for phrase reads: %d sec\n"%(t3 - t2)
    graph = patterns.overlapGraph(k_mer)
    t4 = time.time()
    print "Running time for optimized algorithm: %d sec\n"%(t4 - t2)
    numberOfNodes = len(graph)
    numberOfEdges = sum([len(edges) for edges in graph.values()])
    print "Q3: the total edges are %d\n"%numberOfEdges
    print "Q4: the total nodes are %d"%numberOfNodes
    