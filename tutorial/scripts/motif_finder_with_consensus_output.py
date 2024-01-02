import sys

sys.path.append('C:/Users/samcl/software')

from MotiF.MotiF.motif import *

if __name__ == '__main__':


    mot1 = MotiF('AnnA', output_dir = snakemake.params.dir_out)


    mot1.seq_input(snakemake.input[0])

    mot1.motif_finder()

    for i,consensus in enumerate(mot1.get_motif_consensus(to_fasta=True)):

        print(mot1.get_seq_names(i))
        print(consensus)

    







    