#!/usr/bin/env python
"""
MotiF
==========
    This module has been designed to search DNA and RNA sequences for 
    matches to a motif.
    
"""
import os
import sys
from Bio import SeqIO
import numpy as np


__author__  = 'Sam Clark <sam.clark@york.ac.uk>'
__license__ = 'MIT'
__version__ = '1.0'

class motif_pack():
	"""
    	This class implements 3 main methods:
        load_sequence(self,seq_file): Loads the sequence to search 
        given a file path to a fasta file (seq_file).
        motif_finder(self): Searches the sequence for matches to the 
 	motif provided.
        consensus_generator(self): creates a consensus motif based 
	on the matches the motif provided.
    	"""
        def __init__(self,motif,cl=75):

	    """ Initialises motif_pack object given the motif provided
		and a confidence level for the consensus.
    	    """
            self.motif=motif
            self.consensus=''
            self.cl=cl
            self.match=[]
            self.seq_name=''
            self.seq=''
            self.i=0
            self.ss=''
            self.cA=np.zeros([len(motif)])
            self.cU=np.zeros([len(motif)])
            self.cG=np.zeros([len(motif)])
            self.cC=np.zeros([len(motif)])
            self.cY=np.zeros([len(motif)])
            self.cR=np.zeros([len(motif)])

        def ms_checker(self):

	    """ Checks to see if the notation used for the motif is 
		correct.
    	    """
            library=['A','T','C','G','X','R','U','Y']
            self.motif=self.motif.upper()

            self.motif=self.motif.replace('T','U')
            self.seq=self.seq.transcribe()
                
            if not all(i in library for i in self.motif):
                print self.motif
                print 'Incorrect notation used for motif!'
                quit()
                
        def load_sequence(self,seq_file,select=0):

	    """ Loads the sequence to be searched given a file path
		to a fasta file containing the sequence. 
	    """
            for i,record in enumerate(SeqIO.parse(seq_file,'fasta')):
                if i==select:
                
	            self.seq = record.seq
                    self.seq_name=record.id
            self.ms_checker()

        def match_bin(self):

	    """ Stores match information into numpy array. 
	    """
            if len(self.match)==0:
                self.match=np.array([self.i,str(self.ss)])
            else:
                self.match=np.vstack((self.match,np.array([self.i,str(self.ss)])))
	
	def match_motifs(self):

	    """ Searches sequence for matches to the motif. 
	    """
	    count=0
            motif=self.motif
	    seq=self.ss
            for i in range(len(seq)):
                if motif[i]=='A' and seq[i]=='A':
                    count+=1
                elif motif[i]=='T' and seq[i]=='T':
                    count+=1
                elif motif[i]=='U' and seq[i]=='U':
                    count+=1
                elif motif[i]=='C' and seq[i]=='C':
                    count+=1
                elif motif[i]=='G' and seq[i]=='G':
                    count+=1
                elif motif[i]=='R' and (seq[i]=='G' or seq[i]=='A'):
                    count+=1
                elif motif[i]=='Y' and (seq[i]=='C' or seq[i]=='U'):
                    count+=1
                elif motif[i]=='X':
                    count+=1
            #print count
            if count==len(seq):
                self.match_bin()
        
 	def tally(self):

	    """ Counts the occurences of nucleotides and nucleotide types in
		each position of the motif matches. 
	    """
	    count=0
            motif=self.motif
	    seq=self.ss

            for j in range(len(self.match)):
                select=self.match[j,1]
                for i in range(len(select)):
                    print select
                    if select[i]=='A':
                        self.cA[i]=self.cA[i]+1
                    elif select[i]=='U':
                        self.cU[i]=self.cU[i]+1
                    elif select[i]=='C':
                        self.cC[i]=self.cC[i]+1
                    elif select[i]=='G':
                        self.cG[i]=self.cG[i]+1
                    elif select[i]=='A' or select[i]=='G':
                        self.cR[i]=self.cR[i]+1
                    elif select[i]=='C' or select[i]=='U':
                        self.cY[i]=self.cY[i]+1
        
        def consensus_generator(self):

	    """Creates a consensus sequence based on the motif matches.
	    """
            rat=self.cl/100.0
            self.cA=self.cA/float(len(self.match))
            self.cU=self.cU/float(len(self.match))
            self.cG=self.cG/float(len(self.match))
            self.cC=self.cC/float(len(self.match))
            self.cR=self.cR/float(len(self.match))
            self.cY=self.cY/float(len(self.match))
           
            for i in range(len(self.motif)):
                if self.cA[i]>rat:
                    self.consensus+='A'
                elif self.cU[i]>rat:
                    self.consensus+='U'
                elif self.cG[i]>rat:
                    self.consensus+='G'
                elif self.cC[i]>rat:
                    self.consensus+='C'
                elif self.cR[i]>rat:
                    self.consensus+='R'
                elif self.cY[i]>rat:
                    self.consensus+='Y'
                else:
                    self.consensus+='X'

            print self.consensus

	def motif_finder(self):
	    """Wrapper function for finding matches to the motif provided.
	    """
	
	    if self.motif=='' or self.motif=='':
                print 'sequence or motif not found!'
                exit()
            else:
	    
                window=len(self.motif)
            
                for i in range(len(self.seq)-window):
                    ss=self.seq[i:(i+window)]
                    self.ss=ss
                    self.i=i
                    self.match_motifs()
                    
                            




