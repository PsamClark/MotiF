import sys

sys.path.append('C:/Users/samcl/software')

from MotiF.MotiF.motif import *

if __name__ == '__main__':


    mot1 = MotiF('AnnA',output_dir=snakemake.params.dir_out)


    mot1.seq_input(snakemake.input[0])

    mot1.stats_output()








    