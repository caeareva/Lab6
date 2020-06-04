#!/usr/bin/env python3
# Student: Carlos Arevalo (caeareva)
# Group: None

#####################################################################################################
# File: findUniques.py 
# Executable in stdin: python findUniques.py < [input-file].fa
#                      python findUniques.py < bos-tRNA.fa
#
# Pupose: find unique sequences in tRNAs 
#
#####################################################################################################

'''
In this class, we define objects to read FastaA files inside the sequenceAnalysis program.
'''

import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

''' 
In this class, we find unique sequences only found in a specific tRNA sequence
'''

import sys 
import copy

class tRNA:
    '''
    Class receives an input fasta file, and clean clean all non wanted elements 
    inside the sequences such as dashes, periods and spaces. The method uses FastAreader class 
    to read the input fastas and sets to store the unique elements found in each tRNA. 
    '''
    
    def __init__(self):
        '''
        Initialize program and dictionaries
        '''
        self.powerSetCount = [] # holds power set counts
        self.uniquetRNAs = []  # hold unique power sets
        self.headerSeq = {}  # stores headers and tRNA sequences
        self.rnaSet = set() # # holds tRNA subsets

    def cleanSequence(self, sequence):
        '''
        Remove periods, dashes and underscores from sequence
        '''
        cleanSeq = sequence.replace('.', '').replace('-', '').replace('_', '')
        return cleanSeq

    def getPowerSet(self, sequence):
        '''
        Get power set of sequence.
        Input: tRNA sequence
        Return: a list of subsets form the power sets containing essentials/uniques sequences
        Example:
                    PowerSet("ACG") == {"", "ACG", "A","C", "G", "AC", "CG"}
        '''
        powerSet = set() # power sets for each tRNA
        for subset in range(len(sequence)): # length of tRNA sequence
            lenSeq = len(sequence)
            # return power sets shorter than sequence 
            while lenSeq > subset:
                rna = sequence[subset: lenSeq]
                powerSet.add(rna)  # take tRNA seq and form sets from big to small.
                lenSeq -= 1

        return powerSet

    def findUniques(self):
        '''
        Find unique subsets in each sequence by removing union of all powersets
        '''
        inSequence = FastAreader() # read sequence 
        count = 0  
        for header, sequence in inSequence.readFasta():
            # remove unwated characters in sequence 
            filteredSeq = self.cleanSequence(sequence)
            # builds headerSeq dictionary
            self.headerSeq[count] = [header, filteredSeq]
            mySet = self.getPowerSet(filteredSeq)
            # appends each tRNA sequence to a list
            self.powerSetCount.append(mySet)
            count = count + 1 # same as count += 1
            pass
        
        # get union sets and remove all elements of another set
        for element in self.powerSetCount:
            unionSet = set() # holds union sets
            setList = list(self.powerSetCount) # copy original power sets to a list
            self.rnaSet = element.copy() # copy power sets to list
            setList.remove(self.rnaSet) # removes elements from list
            for item in setList:
                unionSet = (unionSet|item) # perform union of unionSet and item
                # return set power sets after removing elements found in other set
                self.rnaSet.difference_update(unionSet) 
                pass
        
            # compare all unique subsets and find minimum subsets
            essentials = self.rnaSet.copy()
            # remove subsets from tRNAs subsets 
            for subSeq in self.rnaSet:
                uniques = self.rnaSet.copy()
                uniques.remove(subSeq)
                # get essentials in uniques 
                for seq in uniques:
                    if subSeq in seq:
                        # find the minimum subsets
                        if len(subSeq) < len(seq):
                            essentials.discard(seq)
                            
            self.uniquetRNAs.append(essentials)

################################################################################################
# Main program
################################################################################################

def main(inCL=None):
    '''
    Execute program functions, find unique and essential elements in tRNAs
    and build output structure 
    '''
    mytRNA = tRNA()
    mytRNA.findUniques()
    head = mytRNA.headerSeq
    # get and print tRNA header and sequence
    for item in range(0, len(head)):
        rnaInput = head[item]
        header = rnaInput[0] 
        name = header.replace(' ', '')
        print(name)        
        rnaSeq = rnaInput[1] 
        rnaSequence = rnaSeq.replace(' ', '')
        print(rnaSequence ) # print tRNA sequence
        # get and print essentials from each tRNA sequence
        for position in range(0, len(rnaSeq)):
            for element in mytRNA.uniquetRNAs[item]:
                uniqueLength = len(element)
                if element == rnaSeq[position: position + uniqueLength]:
                    uniqueElements = (".")*position + element
                    print(uniqueElements)

if __name__ == "__main__":
    main()

    