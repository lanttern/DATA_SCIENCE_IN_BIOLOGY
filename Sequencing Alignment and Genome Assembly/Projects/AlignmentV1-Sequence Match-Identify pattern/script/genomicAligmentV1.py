# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:42:50 2015
This project includes two classes: findPattern and checkQuality. The findPattern
class finds the occurance and position of a given pattern in a given genomic 
sequence. The checkQuality class exams quality of sequencing for each cycle

@author: zhihuixie
"""


class findPattern ():
    """
    This class finds the occurance and position of a given pattern in a given 
    genomic sequence in a file.
    """
    def __init__(self, pattern, filename):
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
        
    
    def reverseComplement (self):
        """
        generate reverse complement sequence for a given dna sequence
        """
        complement = {"A": "T", "C": "G", "T": "A", "G": "C"}
        revComPattern = "" # reversed compliment pattern
        for nt in self.pattern:
            revComPattern = complement[nt] + revComPattern
        
        return revComPattern
        
    def match(self, string1, string2, numOfMismatch):
        """
        return True or False for matching results of two strings under the offset
        of numOfMismatch
        """
        counter = 0
        if len(string1) != len(string2):
            return False
        # loop over string to compare character
        for i in range(len(string1)):
            if string1[i] != string2[i]:
                counter += 1
        if counter > numOfMismatch:
            return False
        return True
        
    def patternIdentifier (self, numOfMismatch):
        """
        find positions of a given pattern and the reversed complement pattern 
        in a given genome
        """
        patternLength = len(self.pattern)
        genome = self.readGenome()
        revComPattern = self.reverseComplement()
        occurances = []
        
        for i in range (patternLength): # loop over pattern index
            # loop over genome patterns
            for j in range (i, len(genome), patternLength): 
                genomeMotif = genome[j: j+patternLength]
                # compare genomic motif and patterns
                if numOfMismatch == 0:
                    if (self.match(genomeMotif, self.pattern, 0) or \
                        self.match(genomeMotif, revComPattern, 0))\
                       and j not in occurances: # avoid duplicate records
                        occurances.append(j)
                else:
                    if self.match(genomeMotif, self.pattern, numOfMismatch)\
                       and j not in occurances: # avoid duplicate records:
                        occurances.append(j)
        return occurances
    
class checkQuality ():
    """
    The checkQuality class exams quality of sequencing for each cycle
    """
    def __init__ (self, filename):
        self.filename = filename
    
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
        
    def phre33ToQ (self, qualString):
        """
        transform quality string to quality base score
        """
        qScore = []
        for qual in qualString:
            qScore.append(ord(qual) - 33)
        return qScore
    def findPoorQuality(self):
        """
        find the index of poorest q score in each sequencing
        """
        _, qualities = self.readFastq()
        lowestQScoreIndex = []
        for qualString in qualities:
            qScore = self.phre33ToQ(qualString)
            lowestQScoreIndex.append(qScore.index(min(qScore)))
        return lowestQScoreIndex
        
    def countPoorQuality(self):    
        """
        count number of poorest q score in each cycle
        """
        import collections
        return collections.Counter(self.findPoorQuality())
        
    def plotHist(self):
        """
         show the distribution of poorest q score
        """
        import matplotlib.pyplot as plt
        data = self.countPoorQuality()
        plt.bar(data.keys(), data.values())
        plt.show()


if __name__ == "__main__":
    #Test
    filename = "../data/phix.fa"
    pattern = "ATTA"
    patterns = findPattern(pattern,filename)
    print "Test dataset results - Occurances and leftmost offset: "
    print len(patterns.patternIdentifier(0)), min(patterns.patternIdentifier(0)), "\n"
    
    #Questions 1-6
    filename1 = "../data/lambda_virus.fa"
    #Q1: How many times does AGGT or its reverse complement (ACCT) occur in the 
    #lambda virus genome? E.g. if AGGT occurs 10 times and ACCT occurs 12 times, 
    #you should report 22.
    pattern = "AGGT"
    patterns = findPattern(pattern,filename1)
    print "Q1: The 'AGGT' or 'ACCT' occurs %d times \n" \
           %len(patterns.patternIdentifier(0))
    
    #Q2: How many times does TTAA or its reverse complement occur in the lambda 
    #virus genome? Hint: TTAA and its reverse complement are equal, so remember 
    #not to double count.    
    pattern = "TTAA"
    patterns = findPattern(pattern,filename1)
    print "Q2: The 'TTAA' occurs %d times \n" \
           %len(patterns.patternIdentifier(0))
    
    #Q3: What is the offset of the leftmost occurrence of ACTAAGT or its reverse
    #complement in the Lambda virus genome? E.g. if the leftmost occurrence of 
    #ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse 
    #complement ACTTAGT is at offset 29, then report 29.
    pattern = "ACTAAGT"
    patterns = findPattern(pattern,filename1)
    print "Q3: The offset of the leftmost occurrence of ACTAAGT is %d \n" \
           %min(patterns.patternIdentifier(0))
    
    #Q4: What is the offset of the leftmost occurrence of AGTCGA or its reverse 
    #complement in the Lambda virus genome?
    pattern = "AGTCGA"
    patterns = findPattern(pattern,filename1)
    print "Q4: The offset of the leftmost occurrence of AGTCGA is %d \n" \
           %min(patterns.patternIdentifier(0))
           
    #Q5: How many times does TTCAAGCC occur in the Lambda virus genome when 
    #allowing up to 2 mismatches?    
    pattern = "TTCAAGCC"
    patterns = findPattern(pattern,filename1)
    print "Q5: The 'TTCAAGCC' occurs %d times with mismatch less than 2 \n" \
           %len(patterns.patternIdentifier(2))
       
    #Q6: What is the offset of the leftmost occurrence of AGGAGGTT in the 
    #Lambda virus genome when allowing up to 2 mismatches?   
    pattern = "AGGAGGTT"
    patterns = findPattern(pattern,filename1)
    print "Q6: The offset of the leftmost occurrence of AGGAGGTT with mismatch less than 2 is %d \n" \
          %min(patterns.patternIdentifier(2))      
      
    #Q7: Report which sequencing cycle has the problem.    
    filename2 = "../data/ERR037900_1.first1000.fastq"
    qualities = checkQuality(filename2)
    counters = qualities.countPoorQuality()
    print "Q7: The cycle has most frequent poor quality is %d" \
          %max(counters, key = lambda x: counters[x])
    