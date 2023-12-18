import sys

sys.path.append('C:/Users/samcl/software')

from MotiF.MotiF.motif import *

if __name__ == '__main__':


    mot1 = MotiF('AnnA')


    mot1.seq_input(snakemake.input[0])

    mot1.motif_finder()

    for i,matches in enumerate(mot1.get_motif_matches()):

        print(mot1.get_seq_names(i))
        print(matches)







    