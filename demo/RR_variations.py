import os
import sys
from Bio import SeqIO
import numpy as np

if __name__ == '__main__':
   
    import MotiF

    #open output file
    outfile=open('oc43_RR_variations.csv','w')

    #setup motifs
    submotifs=['UAAU','UAGU','UGAU','UGGU']
    
    outfile.write('sub1,sub2,count\n')

    #iterate through submotif combinations
    for i in range(4):
	for j in range(4):
		
		mid=''.join(['x']*8)
		#initialise motif object
		motif = MotiF.motif_pack(submotifs[i]+mid+submotifs[j],75)
		
		#load sequence    
		motif.load_sequence('HCoV_OC43_ref.fa')
		
		#find motif matches
		motif.motif_finder()
	
		#write out matchs to outfile
		if isinstance(motif.match,list):
		  
			print i,0
			print motif.match
    			outfile.write('%s,%s,%d\n'%(submotifs[i],submotifs[j],0))

		elif len(list(motif.match.shape))==1:
			print i,1
			print motif.match
    			outfile.write('%s,%s,%d\n'%(submotifs[i],submotifs[j],1))
		else:	
			print i,len(motif.match)
			print motif.match
    			outfile.write('%s,%s,%d\n'%(submotifs[i],submotifs[j],len(motif.match)))

    outfile.close()
