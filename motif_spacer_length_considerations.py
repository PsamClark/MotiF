import os
import sys
from Bio import SeqIO
import numpy as np

if __name__ == '__main__':
   
    import MotiF

    #open output file
    outfile=open('oc43_spacer_control.csv','w')
 
    outfile.write('spacer,count\n')
    
    #iterate through spacers of different lengths
    for i in range(41):
	#set spacer
	mid=''.join(['x']*i)
	
	#initialise motif_pack object
	motif = MotiF.motif_pack('RUUR'+mid+'RUUR',75)
	
	#load sequence    
	motif.load_sequence('HCoV_OC43_ref.fa')

	#search for motif matches
	motif.motif_finder()
	  
	#output number of motif matches for each spacer length. 
    	outfile.write('%d,%d\n'%(i,len(motif.match)))

    outfile.close()