#!/usr/bin/env python
"""
MotiF
==========
    This module has been designed to search DNA and RNA sequences for
    matches to a motif.
"""
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd

__author__ = 'Sam Clark <sam.clark@york.ac.uk>'
__license__ = 'MIT'
__version__ = '2.0'

class MotifError(Exception):
    pass

class InputError(Exception):
    pass

class FileError(Exception):

    pass

class SeqError(Exception):
    pass

class MotiF:
    """Motif detection class

    Attributes:
    motif (str): Motif to search for in the sequences.
    _cl (float): Cutoff limit for consensus generation.
    _sig (bool): Logical indicator of whether significance is being
    calculated.
    match (:obj: 'list' of :obj: 'np.ndarray'): List of matches to
    search motif.
    seq_names (:obj: 'list' of :obj: 'str'): Names of sequences under
    investigation.
    seqs (:obj: 'list' of :obj: 'str'): Sequences under investigation.
    consensus (str): Consensus motif for matches.
    _cA (np.array): Frequency of A in motif matches.
    _cG (np.array): Frequency of G in motif matches.
    _cC (np.array): Frequency of C in motif matches.
    _cT (np.array): Frequency of T in motif matches.
    _cU (np.array): Frequency of U in motif matches.
    _cY (np.array): Frequency of Y in motif matches.
    _cR (np.array): Frequency of R in motif matches.
    _cM (np.array): Frequency of M in motif matches.
    _cK (np.array): Frequency of K in motif matches.
    _cS (np.array): Frequency of S in motif matches.
    _cW (np.array): Frequency of W in motif matches.
    _cH (np.array): Frequency of H in motif matches.
    _cB (np.array): Frequency of B in motif matches.
    _cV (np.array): Frequency of V in motif matches.
    _cD (np.array): Frequency of D in motif matches.
    _aval (:obj: 'list' of :obj: 'int'): Number of times motif matched in
    target sequence.
    _raw_vals (:obj: 'list' of :obj: 'list' of :obj: 'int'):Number of
    times motif matched in each of randomised sequences.
    _std (:obj: 'list' of :obj: 'float'): Standard deviation of occurences
    of motif in randomised sequences.
    _mean (:obj: 'list' of :obj: 'float'): Mean of occurences of motif in
    randomised sequences.
    zscore: (:obj: 'list' of :obj: 'float'): Z-score for motif count in
    sequence under investigation.
    _rand_seqs (:obj: 'list' of :obj: 'list' of :obj: 'int'): Randomised
    sequences.
    _acont (:obj: 'list' of :obj: 'float'): A content of sequences.
    _gcont (:obj: 'list' of :obj: 'float'): G content of sequences.
    _ccont (:obj: 'list' of :obj: 'float'): C content of sequences.
    _tcont (:obj: 'list' of :obj: 'float'): T content of sequences.
    _ucont (:obj: 'list' of :obj: 'float'): U content of sequences.
    _ycont (:obj: 'list' of :obj: 'float'): Y content of sequences.
    _rcont (:obj: 'list' of :obj: 'float'): R content of sequences.
    _mcont (:obj: 'list' of :obj: 'float'): M content of sequences.
    _kcont (:obj: 'list' of :obj: 'float'): K content of sequences.
    _scont (:obj: 'list' of :obj: 'float'): S content of sequences.
    _wcont (:obj: 'list' of :obj: 'float'): W content of sequences.
    _hcont (:obj: 'list' of :obj: 'float'): H content of sequences.
    _bcont (:obj: 'list' of :obj: 'float'): B content of sequences.
    _vcont (:obj: 'list' of :obj: 'float'): V content of sequences.
    _dcont (:obj: 'list' of :obj: 'float'): D content of sequences.
    _ncont (:obj: 'list' of :obj: 'float'): N content of sequences.
   """

    def __init__(self, motif='', cl=75, output_dir = None):
        """ Initials MotiF class.
        Args:
        motif (str): Motif to search for in sequences under investigation.
        cl (float): Cutoff percentage for consensus generation.
        """

        # load and check motif

        if not isinstance(motif,str):
            raise MotifError('motif must be a string!')

        elif len(motif) == 0:
            raise MotifError('motif must be supplied!')
            
        self._motif = motif

        self._motif_checker()

        # initialise and load parameters

        if (not isinstance(cl,int) and not isinstance(cl,float)) or \
            (isinstance(cl,bool) or  not np.isfinite(cl)):
            raise InputError('consensus limit must be float or integer!')

        self._cl = cl
        self._sig = False
        self._seq_file = ''

        # initialise seqeunce and motif variables
        self._matches = []
        self._seq_names = []
        self._seqs = []
        self._consensus = []
        self._dir_out = ''

        #if output_dir specified create it
        if output_dir is not None:
            self.set_dir_out(output_dir)

        # initialise temporary storage variables
        self._i = 0
        self._ss = ''
        self._temp_match = []

        # initialise motif stats
        self._cA = []
        self._cU = []
        self._cT = []
        self._cG = []
        self._cC = []
        self._cY = []
        self._cR = []
        self._cM = []
        self._cK = []
        self._cS = []
        self._cW = []
        self._cH = []
        self._cB = []
        self._cV = []
        self._cD = []

        # initialise significance variables
        self._aval = []
        self._raw_vals = []
        self._ps_mval = []
        self._ps_mval_rand = []
        self._ps_mean = []
        self._ps_std = []
        self._ps_zscore = []
        self._nruns = 100
        self._std = []
        self._zscore = []
        self._mean = []
        self._rand_seqs = []

        # initialise sequence statistic variables
        self._acont = []
        self._ccont = []
        self._gcont = []
        self._tcont = []
        self._ucont = []

        self._rcont = []
        self._ycont = []
        self._ncont = []
        self._mcont = []
        self._kcont = []
        self._scont = []
        self._wcont = []
        self._hcont = []
        self._bcont = []
        self._vcont = []
        self._dcont = []

    def set_dir_out(self,new_path,delete_previous = False):
        """Output data into new directory. 

        Args:

            new_path (str): Path of new directory.
            delete_previous (bool): Request if previous output directory
                should be deleted.
        """

        #delete previous directory if requested. 
        if delete_previous:
            shutil.rmtree(self._dir_out)
        
        #set new dir out path.
        self._dir_out = new_path

        if new_path[-1] != '/':
            self._dir_out += '/'

        #create dir out if it doesn't exist.
        if not os.path.exists(self._dir_out):
            os.mkdir(self._dir_out)

    def get_dir_out(self):
        """retrieve dir out path
        """
        #print dir out path
        print(self._dir_out)

    def _motif_checker(self):
        """ Checks to see if the notation used for the motif is
            correct.
        """

        # specify the library of the IUPAC-IUB code
        library = ['A', 'T', 'C', 'G', 'R', 'U', 'Y', 'N', 'M',
                   'K', 'S', 'W', 'H', 'B', 'V', 'D']

        # convert characters in motif to uppercase
        self._motif = self._motif.upper()

        # ensure that motif complies with IUPAC-IUB code
        if not all(i in library for i in self._motif):
            raise MotifError('Please use IUPAC-IUB code for motif!')
            

    def seq_input(self, seq_file=''):
        """ Loads the sequence to be searched given a file path
            to a fasta file containing the sequence.
        Args:
            seq_file (str): Sequence FASTA file name/path.
        """

        if not isinstance(seq_file,str):
            raise FileError('Sequence file name/path must be string!')

        elif len(seq_file) == 0:
            raise FileError('Sequence file must be provided!')

        # ierate through sequence file and extract sequence information
        
        if '.fastq' in seq_file: 
            with open(seq_file,'r') as handle:
                seq_recs=SeqIO.parse(handle, "fastq")
                if not any(seq_recs):
                    raise FileError('File not in FASTQ format!')
                
            
            for record in SeqIO.parse(seq_file, 'fastq'):
                seq=str(record.seq).upper()                
                if len(seq) == 0:
                    raise SeqError('Empty sequence detected in record: ' +
                                record.id)
                elif seq.count('T')>0 and seq.count('U')>0:
                    raise SeqError('Sequence ' + record.id +
                                ' contains both U and T bases!')
                self._seq_names.append(str(record.id).replace(' ', ''))
                self._seqs.append(str(record.seq).upper())
            if len(self._seqs) == 0:
                raise FileError('File given had no sequences in it!')
        
        else: 
                
            with open(seq_file,'r') as handle:
                seq_recs=SeqIO.parse(handle, "fasta")
                if not any(seq_recs):
                    raise FileError('File not in FASTA format!')

            for record in SeqIO.parse(seq_file, 'fasta'):
                seq=str(record.seq).upper()                
                if len(seq) == 0:
                    raise SeqError('Empty sequence detected in record: ' +
                                record.id)
                elif seq.count('T')>0 and seq.count('U')>0:
                    raise SeqError('Sequence ' + record.id +
                                ' contains both U and T bases!')
                self._seq_names.append(str(record.id).replace(' ', ''))
                self._seqs.append(str(record.seq).upper())
            if len(self._seqs) == 0:
                raise FileError('File given had no sequences in it!')
            
            # convert ambiguous characters in sequences to N
            self._ambig_decider()

            # calculate statistics for sequence
            self._seq_stats()
            
            self._set_tally_vars()
        
        self._seq_file = seq_file

    def _set_tally_vars(self):
        """ Set motif stat variables
        """

        # initialise motif stats
        self._cA = np.zeros([len(self._motif),len(self._seqs)])
        self._cU = np.zeros([len(self._motif),len(self._seqs)])
        self._cT = np.zeros([len(self._motif),len(self._seqs)])
        self._cG = np.zeros([len(self._motif),len(self._seqs)])
        self._cC = np.zeros([len(self._motif),len(self._seqs)])
        self._cY = np.zeros([len(self._motif),len(self._seqs)])
        self._cR = np.zeros([len(self._motif),len(self._seqs)])
        self._cM = np.zeros([len(self._motif),len(self._seqs)])
        self._cK = np.zeros([len(self._motif),len(self._seqs)])
        self._cS = np.zeros([len(self._motif),len(self._seqs)])
        self._cW = np.zeros([len(self._motif),len(self._seqs)])
        self._cH = np.zeros([len(self._motif),len(self._seqs)])
        self._cB = np.zeros([len(self._motif),len(self._seqs)])
        self._cV = np.zeros([len(self._motif),len(self._seqs)])
        self._cD = np.zeros([len(self._motif),len(self._seqs)])




    def get_seq_names(self,i=None):

        """Get the names of sequences
        """
        
        #if specific sequence name requested return that one. 
        if i is not None: 

            return self._seq_names[i]
        
        #else return all names.     
        else:

            return self._seq_names

    def get_motif_matches(self):
        """ get the matches for the motifs in the different sequences.
        """
        return self._matches
    

    def _match_bin(self):
        """ Stores match information into numpy array.
        """
        # if statement based on whether significances are being calculated
        if self._sig:

            # if nothing in temp_match storage initialise
            if len(self._temp_match) == 0:
                self._temp_match = pd.DataFrame([[self._i+1, str(self._ss)]],
                        columns=['position','motif'],index=['0'])

            # else add to preexisting temp match object
            else:
                index = int(self._temp_match.index[-1])+1
                self._temp_match.loc[str(index)] = \
                        [self._i+1, str(self._ss)]
        # else perform same task using the match variable in self
        else:

            if len(self._matches[-1]) == 0:
                self._matches[-1] = pd.DataFrame([[self._i+1, str(self._ss)]],
                        columns=['position','motif'],index=['0'])
            else:

                index = int(self._matches[-1].index[-1])+1

                self._matches[-1].loc[str(index)] = [self._i+1,
                        str(self._ss)]

    def _match_motifs(self):
        """ Searches sequence for matches to the motif.
        """
        # initialise counter
        count = 0
        # extract motif and sub sequence to compare
        motif = self._motif
        seq = self._ss

        # iteratively score the subsequence based on how each character matchs
        # the corresponding character in the motif
        for i in range(len(seq)):
            if motif[i] == 'A' and seq[i] == 'A':
                count += 1
            elif motif[i] == 'T' and seq[i] == 'T':
                count += 1
            elif motif[i] == 'U' and seq[i] == 'U':
                count += 1
            elif motif[i] == 'C' and seq[i] == 'C':
                count += 1
            elif motif[i] == 'G' and seq[i] == 'G':
                count += 1
            elif motif[i] == 'R' and (seq[i] == 'G' or seq[i] == 'A'):
                count += 1
            elif motif[i] == 'Y' and (seq[i] == 'C' or (seq[i] == 'U' or
                                                        seq[i] == 'T')):
                count += 1
            elif motif[i] == 'M' and (seq[i] == 'A' or seq[i] == 'C'):
                count += 1
            elif motif[i] == 'K' and (seq[i] == 'G' or (seq[i] == 'U' or
                                                        seq[i] == 'T')):
                count += 1
            elif motif[i] == 'S' and (seq[i] == 'G' or seq[i] == 'C'):
                count += 1
            elif motif[i] == 'W' and (seq[i] == 'A' or
                                      (seq[i] == 'U' or seq[i] == 'T')):
                count += 1
            elif motif[i] == 'H' and ((seq[i] == 'A' or seq[i] == 'C') or
                                      (seq[i] == 'U' or seq[i] == 'T')):
                count += 1
            elif motif[i] == 'B' and ((seq[i] == 'G' or seq[i] == 'C') or
                                      (seq[i] == 'U' or seq[i] == 'T')):
                count += 1
            elif motif[i] == 'V' and (seq[i] == 'A' or
                                      (seq[i] == 'C' or seq[i] == 'G')):
                count += 1
            elif motif[i] == 'D' and ((seq[i] == 'A' or seq[i] == 'G') or
                                      (seq[i] == 'U' or seq[i] == 'T')):
                count += 1
            elif motif[i] == 'N':
                count += 1

        # if count is the same length as the motif bin the match
        if count == len(seq):
            self._match_bin()

    def _tally(self):
        """ Counts the occurences of nucleotides and nucleotide types in
            each position of the motif matches.
        """
        print(self._cR)

        # iterate through each set of matches
        for k in range(len(self._matches)):

            # iterate through the matches
            for j in range(len(self._matches[k])):
                select = self._matches[k].loc[str(j),'motif']

                # calculate contents of the matches at each position
                for i in range(len(select)):
                    if select[i] == 'A':
                        self._cA[i,k] += 1
                    elif select[i] == 'U':
                        self._cU[i,k] += 1

                    elif select[i] == 'T':
                        self._cT[i,k] += 1
                    elif select[i] == 'C':
                        self._cC[i,k] += 1
                    elif select[i] == 'G':
                        self._cG[i,k] += 1
                    if select[i] == 'A' or select[i] == 'G':
                        self._cR[i,k] += 1
                    if select[i] == 'C' or (select[i] == 'U' or
                                              select[i] == 'T'):
                        self._cY[i,k] += 1
                    if select[i] == 'A' or select[i] == 'C':
                        self._cM[i,k] += 1
                    if select[i] == 'G' or (select[i] == 'U' or
                                              select[i] == 'T'):
                        self._cK[i,k] += 1
                    if select[i] == 'C' or select[i] == 'G':
                        self._cS[i,k] += 1
                    if select[i] == 'A' or (select[i] == 'U' or
                                              select[i] == 'T'):
                        self._cW[i,k] += 1
                    if (select[i] == 'A' or select[i] == 'C') or \
                    (select[i] == 'U' or select[i] == 'T'):
                        self._cH[i,k] += 1
                    if (select[i] == 'C' or select == 'G') or \
                    (select[i] == 'U' or select[i] == 'T'):
                        self._cB[i,k] += 1
                    if select[i] == 'C' or (select[i] == 'G' or
                                              select[i] == 'A'):
                        self._cV[i,k] += 1
                    if (select[i] == 'G' or select[i] == 'A') or \
                    (select[i] == 'U' or select[i] == 'T'):
                        self._cD[i,k] += 1

            print (len(self._matches[k]))
        print (self._cR)


    def get_motif_consensus(self, to_fasta = False):
        """ Get motif consensuses for the different sequences. 
        """
        #tally the nucleotides in each motif position
        self._tally()

        #generate consensuses. 

        self._consensus_generator()

        if to_fasta:
            self.output_consenses()

        return self._consensus

    def _consensus_generator(self):
        """Creates a consensus sequence based on the motif matches.
        """

        # convert consensus cutoff into ratio
        rat = self._cl / 100.0

        # iterate through the investigated sequences
        for j in range(len(self._seqs)):
            # convert motif stats from counts to percentages
            self._cA[:,j] = self._cA[:,j] / float(len(self._matches[j]))
            self._cU[:, j] = self._cU[:,j] / float(len(self._matches[j]))
            self._cT[:, j] = self._cT[:, j] / float(len(self._matches[j]))
            self._cG[:, j] = self._cG[:, j] / float(len(self._matches[j]))
            self._cC[:, j] = self._cC[:, j] / float(len(self._matches[j]))
            self._cR[:, j] = self._cR[:, j] / float(len(self._matches[j]))
            self._cY[:, j] = self._cY[:, j] / float(len(self._matches[j]))
            self._cM[:, j] = self._cM[:, j] / float(len(self._matches[j]))
            self._cK[:, j] = self._cK[:, j] / float(len(self._matches[j]))
            self._cS[:, j] = self._cS[:, j] / float(len(self._matches[j]))
            self._cW[:, j] = self._cW[:, j] / float(len(self._matches[j]))
            self._cH[:, j] = self._cH[:, j] / float(len(self._matches[j]))
            self._cB[:, j] = self._cB[:, j] / float(len(self._matches[j]))
            self._cV[:, j] = self._cV[:, j] / float(len(self._matches[j]))
            self._cD[:, j] = self._cD[:, j] / float(len(self._matches[j]))
            
            # initialise motif consensus
            consensus = ''

            # iterate across motif selecting characters for the consensus
            # motif based on whether they meet the criteria
            for i in range(len(self._motif)):
                if self._cA[i, j] > rat:
                    consensus += 'A'
                elif self._cU[i, j] > rat:
                    consensus += 'U'
                elif self._cT[i, j] > rat:
                    consensus += 'T'
                elif self._cG[i, j] > rat:
                    consensus += 'G'
                elif self._cC[i, j] > rat:
                    consensus += 'C'
                elif self._cR[i, j] > rat:
                    consensus += 'R'
                elif self._cY[i, j] > rat:
                    consensus += 'Y'
                elif self._cM[i, j] > rat:
                    consensus += 'M'
                elif self._cK[i, j] > rat:
                    consensus += 'K'
                elif self._cS[i, j] > rat:
                    consensus += 'S'
                elif self._cW[i, j] > rat:
                    consensus += 'W'
                elif self._cH[i, j] > rat:
                    consensus += 'H'
                elif self._cB[i, j] > rat:
                    consensus += 'B'
                elif self._cV[i, j] > rat:
                    consensus += 'V'
                elif self._cD[i, j] > rat:
                    consensus += 'D'
                else:
                    consensus += 'N'

            # store the consensus
            self._consensus.append(consensus)

    def output_consenses(self):
        """Output the consenses for the motif in each seq
        """
        
        with open(self._dir_out + self._motif+'_consenses.fa', 'w') as file_out:

            for i in range(len(self._consensus)):

                file_out.write('>' + self._seq_names[i] + '\n')
                file_out.write(self._consensus[i] + '\n')


    def motif_finder(self, i=None, j=None,
         to_csv = False):
        """Wrapper function for finding matches to the motif provided.
        """
        # if significance calculations being performed do so over the
        # random sequences
        if self._sig:
            self._temp_match = []
            self._mscan(self._rand_seqs[i][j])


        # else perform the scan over the investigated sequences
        else:
            for j in range(len(self._seqs)):
                self._matches.append([])
                self._mscan(self._seqs[j])

            # output matches to file
            if to_csv:
                self._data_output()

    def _mscan(self, seq):
        """Scan through the sequence for motif.
        Args:
           seq (str): sequence to scan across
        """

        # set the window length
        window = len(self._motif)

        # window along the sequence
        for i in range(len(seq) - window + 1):
            ss = seq[i:(i + window)]
            self._ss = ss
            self._i = i
            self._match_motifs()
            
    def _seq_stats(self):
        """Calculate content statistics of sequences.
        """
        # iterate through investigated sequences
        for i in self._seqs:
            # count the number of nucleotides in sequences
            self._acont.append(i.count('A') / float(len(i)))
            self._gcont.append(i.count('G') / float(len(i)))
            self._ccont.append(i.count('C') / float(len(i)))
            self._tcont.append(i.count('T') / float(len(i)))
            self._ucont.append(i.count('U') / float(len(i)))
            self._ncont.append(i.count('N') / float(len(i)))
            
            # calculate composite statistics
            self._rcont.append(self._acont[-1] + self._gcont[-1])
            self._ycont.append(1 - self._rcont[-1])
            self._mcont.append(self._acont[-1] + self._ccont[-1])
            self._kcont.append(
                self._gcont[-1] + self._ucont[-1] + self._tcont[-1])
            self._scont.append(self._ccont[-1] + self._gcont[-1])
            self._wcont.append(self._acont[-1] + self._tcont[-1] + \
                self._ucont[-1])
            self._hcont.append(
                self._acont[-1] + self._ccont[-1] + \
                self._ucont[-1] + self._tcont[-1])
            self._bcont.append(
                self._gcont[-1] + self._ccont[-1] + \
                self._ucont[-1] + self._tcont[-1])
            self._vcont.append(
                self._acont[-1] + self._ccont[-1] + self._gcont[-1])
            self._dcont.append(
                self._acont[-1] + self._gcont[-1] + \
                self._ucont[-1] + self._tcont[-1])

    def _ambig_decider(self):
        """Convert ambiguous nucleotides in sequences to N.
        """

        # accepted nucleotides
        nucs = ['A', 'C', 'G', 'T', 'U', 'N']

        # iterate through sequences
        for i in range(len(self._seqs)):

            # iterate along the sequence
            for j in range(len(self._seqs[i])):

                # if nucleotide is not acceptable convert to N
                if self._seqs[i][j] not in nucs:
                    self._seqs[i] = self._seqs[i][:j] + \
                        'N' + self._seqs[i][j + 1:]

    def _seq_randomiser(self):
        """Generate randomised sequences based on composition of investigated
        sequences.
        """

        # acceptable nucleotides
        nucs = ['A', 'C', 'G', 'T', 'U', 'N']

        # iterate through the sequences
        for i in range(len(self._seqs)):

            # populate probability array:
            probs = [self._acont[i], self._ccont[i], self._gcont[i],
                     self._tcont[i], self._ucont[i], self._ncont[i]]

            # initiate random sequence dumping array
            self._rand_seqs.append([])

            # generate random sequences
            for j in range(self._nruns):
                self._rand_seqs[-1].append(''.join(np.random.choice(nucs,
                                           len(self._seqs[i]), p=probs)))

    def motif_signif_calculator(self, n=100):
        """Significance calculations wrapper.
        Args:
            n (int): Number of randomised sequences to calculate
            significance over.
        """

        # set parameters
        if not isinstance(n,int) or isinstance(n,bool):
            raise InputError('n must be integer!')
        self._nruns = n
        self._sig = True

        self._seq_randomiser()

        # run through calculations
        self._sig_calculator()

        #output data
        self._sig_output()


    def _sig_calculator(self):
        """Calculate the significance of motif occurence
        """

        # iterate through sequences
        for i in range(len(self._seqs)):

            # set actual value
            if self._matches[i].ndim == 1:
                self._aval.append(1)
            else:
                self._aval.append(len(self._matches[i]))

            # create raw val bin
            self._raw_vals.append([])

            # iterate through random sequences
            for j in range(len(self._rand_seqs[i])):

                # run motif finding algorithm
                self.motif_finder(i, j)

                # dump raw motif match count
                self._raw_vals[-1].append(len(self._temp_match))

            # calculate mean and standard deviations of motif match
            # numbers from randomised sequences
            self._mean.append(np.mean(self._raw_vals[-1]))
            self._std.append(np.std(self._raw_vals[-1]))

            # calculate Z-score of match count on ivestigated sequence
            if self._std[-1] != 0:

                self._zscore.append(
                    (self._aval[-1] - self._mean[-1]) / self._std[-1])
            else:
                self._zscore.append(0)

        # output significance values
        self._sig_output()

    def _sig_output(self):
        """Writes significance values to file

        """
        file_name=self._dir_out+self._motif + '_sig.csv'

        with open(file_name, 'w') as file_out:

            # write headers
            file_out.write('seq,count,mean,std,zscore\n')

            # iterate through significance values and write to file
            for i in range(len(self._seqs)):
                file_out.write(
                    '%s,%d,%f,%f,%f\n' %
                    (self._seq_names[i],
                        self._aval[i],
                        self._mean[i],
                        self._std[i],
                        self._zscore[i]))


    def _data_output(self):
        """Output motif matches to file.
        """

        # iterate through sequences
        for j in range(len(self._seq_names)):

            # write matches to file
            if len(self._matches[j])>0:
                self._matches[j].to_csv(self._dir_out+self._seq_names[j] + '_' + \
                        self._motif + '_occurences.csv', index = False)

    def stats_output(self):
        """output sequence stats to output directory"""

        for j in range(len(self._seq_names)):
            with open(self._dir_out+self._seq_names[j] + '_statistics.csv', 
                      'w') as file_out:
                file_out.write('metric,value (%)\n')
                file_out.write('length,%d\n' % (len(self._seqs[j])))
                file_out.write('GC,%f\n'%((self._gcont[j]+self._ccont[j])*100))
                file_out.write('A,%f\n' % (self._acont[j] * 100))
                file_out.write('G,%f\n' % (self._gcont[j] * 100))
                file_out.write('C,%f\n' % (self._ccont[j] * 100))
                file_out.write('T,%f\n' % (self._tcont[j] * 100))
                file_out.write('U,%f\n' % (self._ucont[j] * 100))
                file_out.write('N,%f\n' % (self._ncont[j] * 100))
                file_out.write('M,%f\n' % (self._mcont[j] * 100))
                file_out.write('K,%f\n' % (self._kcont[j] * 100))
                file_out.write('S,%f\n' % (self._scont[j] * 100))
                file_out.write('W,%f\n' % (self._wcont[j] * 100))
                file_out.write('H,%f\n' % (self._hcont[j] * 100))
                file_out.write('B,%f\n' % (self._bcont[j] * 100))
                file_out.write('V,%f\n' % (self._vcont[j] * 100))
                file_out.write('D,%f\n' % (self._dcont[j] * 100))

