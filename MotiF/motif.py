#!/usr/bin/env python
"""
MotiF
==========
    This module has been designed to search DNA and RNA sequences for
    matches to a motif.
"""
import os
import sys
print(sys.path)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from Twrap.Twrap_v2 import *


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

    def __init__(self, motif='', cl=75):
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
        self.match = []
        self.seq_names = []
        self.seqs = []
        self.consensus = ''

        # initialise temporary storage variables
        self._i = 0
        self._ss = ''
        self._temp_match = []

        # initialise motif stats
        self._cA = np.zeros([len(motif)])
        self._cU = np.zeros([len(motif)])
        self._cT = np.zeros([len(motif)])
        self._cG = np.zeros([len(motif)])
        self._cC = np.zeros([len(motif)])
        self._cY = np.zeros([len(motif)])
        self._cR = np.zeros([len(motif)])
        self._cM = np.zeros([len(motif)])
        self._cK = np.zeros([len(motif)])
        self._cS = np.zeros([len(motif)])
        self._cW = np.zeros([len(motif)])
        self._cH = np.zeros([len(motif)])
        self._cB = np.zeros([len(motif)])
        self._cV = np.zeros([len(motif)])
        self._cD = np.zeros([len(motif)])

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
            self.seq_names.append(str(record.id).replace(' ', ''))
            self.seqs.append(str(record.seq).upper())
        if len(self.seqs) == 0:
            raise FileError('File given had no sequences in it!')
        
        self._seq_file = seq_file
        # convert ambiguous characters in sequences to N
        self._ambig_decider()

        # calculate statistics for sequence
        self._seq_stats()


    def _match_bin(self):
        """ Stores match information into numpy array.
        """
        # if statement based on whether significances are being calculated
        if self._sig:

            # if nothing in temp_match storage initialise
            if len(self._temp_match) == 0:
                self._temp_match = np.array((self._i+1, str(self._ss)))

            # else add to preexisting temp match object
            else:
                self._temp_match = np.vstack(
                    (self._temp_match, np.array([self._i+1, str(self._ss)])))
        # else perform same task using the match variable in self
        else:

            if len(self.match[-1]) == 0:
                self.match[-1] = np.array((self._i+1, str(self._ss)))
            else:
                self.match[-1] = np.vstack((self.match[-1],
                                            np.array([self._i+1,
                                                     str(self._ss)])))

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

        # iterate through each set of matches
        for k in range(len(self.match)):

            # iterate through the matches
            for j in range(len(self.match[k])):
                select = self.match[k][j, 1]

                # calculate contents of the matches at each position
                for i in range(len(select)):
                    print(select)
                    if select[i] == 'A':
                        self.cA[k][i] += 1
                    elif select[i] == 'U':
                        self.cU[k][i] += 1
                    elif select[i] == 'C':
                        self.cC[k][i] += 1
                    elif select[i] == 'G':
                        self.cG[k][i] += 1
                    elif select[i] == 'A' or select[i] == 'G':
                        self.cR[k][i] += 1
                    elif select[i] == 'C' or (select[i] == 'U' or
                                              select[i] == 'T'):
                        self.cY[k][i] += 1
                    elif select[i] == 'A' or select[i] == 'C':
                        self.cM[k][i] += 1
                    elif select[i] == 'G' or (select[i] == 'U' or
                                              select[i] == 'T'):
                        self.cK[k][i] += 1
                    elif select[i] == 'C' or select[i] == 'G':
                        self.cS[k][i] += 1
                    elif select[i] == 'A' or (select[i] == 'U' or
                                              select[i] == 'T'):
                        self.cW[k][i] += 1
                    elif (select[i] == 'A' or select[i] == 'C') or \
                    (select[i] == 'U' or select[i] == 'T'):
                        self.cH[k][i] += 1
                    elif (select[i] == 'C' or select == 'G') or \
                    (select[i] == 'U' or select[i] == 'T'):
                        self.cB[k][i] += 1
                    elif select[i] == 'C' or (select[i] == 'G' or
                                              select[i] == 'A'):
                        self.cV[k][i] += 1
                    elif (select[i] == 'G' or select[i] == 'A') or \
                    (select[i] == 'U' or select[i] == 'T'):
                        self.cD[k][i] += 1

    def _consensus_generator(self):
        """Creates a consensus sequence based on the motif matches.
        """

        # convert consensus cutoff into ratio
        rat = self.cl / 100.0

        # iterate through the investigated sequences
        for j in range(len(self.seqs)):
            # convert motif stats from counts to percentages
            self.cA[j] = self.cA[j] / float(len(self.match[j]))
            self.cU[j] = self.cU[j] / float(len(self.match[j]))
            self.cT[j] = self.cT[j] / float(len(self.match[j]))
            self.cG[j] = self.cG[j] / float(len(self.match[j]))
            self.cC[j] = self.cC[j] / float(len(self.match[j]))
            self.cR[j] = self.cR[j] / float(len(self.match[j]))
            self.cY[j] = self.cY[j] / float(len(self.match[j]))
            self.cM[j] = self.cM[j] / float(len(self.match[j]))
            self.cK[j] = self.cK[j] / float(len(self.match[j]))
            self.cS[j] = self.cS[j] / float(len(self.match[j]))
            self.cW[j] = self.cW[j] / float(len(self.match[j]))
            self.cH[j] = self.cH[j] / float(len(self.match[j]))
            self.cB[j] = self.cB[j] / float(len(self.match[j]))
            self.cV[j] = self.cV[j] / float(len(self.match[j]))
            self.cD[j] = self.cD[j] / float(len(self.match[j]))

            # initialise motif consensus
            consensus = ''

            # iterate across motif selecting characters for the consensus
            # motif based on whether they meet the criteria
            for i in range(len(self._motif)):
                if self.cA[j][i] > rat:
                    consensus += 'A'
                elif self.cU[j][i] > rat:
                    consensus += 'U'
                elif self.cT[j][i] > rat:
                    consensus += 'T'
                elif self.cG[j][i] > rat:
                    consensus += 'G'
                elif self.cC[j][i] > rat:
                    consensus += 'C'
                elif self.cR[j][i] > rat:
                    consensus += 'R'
                elif self.cY[j][i] > rat:
                    consensus += 'Y'
                elif self.cM[j][i] > rat:
                    consensus += 'M'
                elif self.cK[j][i] > rat:
                    consensus += 'K'
                elif self.cS[j][i] > rat:
                    consensus += 'S'
                elif self.cW[j][i] > rat:
                    consensus += 'W'
                elif self.cH[j][i] > rat:
                    consensus += 'H'
                elif self.cB[j][i] > rat:
                    consensus += 'B'
                elif self.cV[j][i] > rat:
                    consensus += 'V'
                elif self.cD[j][i] > rat:
                    consensus += 'D'
                else:
                    consensus += 'N'

            # store the consensus
            self.consensus.append(consensus)

    def motif_finder(self, i=None, j=None,method='ps',
            fold_params={'win_size':None,'win_inc':None,
                         'sl_method':'loops'}):
        """Wrapper function for finding matches to the motif provided.
        """
        # if significance calculations being performed do so over the
        # random sequences
        if self._sig:
            self._temp_match = []
            self._mscan(self._rand_seqs[i][j])

        # else perform the scan over the investigated sequences
        else:
            methods=['linear','ps']

            if not isinstance(method,str):
                raise InputError('method must be string')
            elif method not in methods:
                raise InputError('method not found')

            if method == 'ps':
                seq_rec=[]
                for j in range(len(self.seqs)):
                   seq_rec.append(SeqRecord(Seq(self.seqs[j]),
                       id = self.seq_names[j]))
                self._fold(seq_rec,fold_params)
                self._sl_motif_filter()
            else:
                for j in range(len(self.seqs)):
                    self.match.append([])
                    self._mscan(self.seqs[j])
                # output matches to file
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
        for i in self.seqs:
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

        # output statistics
        self.stats_output()

    def _ambig_decider(self):
        """Convert ambiguous nucleotides in sequences to N.
        """

        # accepted nucleotides
        nucs = ['A', 'C', 'G', 'T', 'U', 'N']

        # iterate through sequences
        for i in range(len(self.seqs)):

            # iterate along the sequence
            for j in range(len(self.seqs[i])):

                # if nucleotide is not acceptable convert to N
                if self.seqs[i][j] not in nucs:
                    self.seqs[i] = self.seqs[i][:j] + \
                        'N' + self.seqs[i][j + 1:]

    def _seq_randomiser(self):
        """Generate randomised sequences based on composition of investigated
        sequences.
        """

        # acceptable nucleotides
        nucs = ['A', 'C', 'G', 'T', 'U', 'N']

        # iterate through the sequences
        for i in range(len(self.seqs)):

            # populate probability array:
            probs = [self._acont[i], self._ccont[i], self._gcont[i],
                     self._tcont[i], self._ucont[i], self._ncont[i]]

            # initiate random sequence dumping array
            self._rand_seqs.append([])

            # generate random sequences
            for j in range(self._nruns):
                self._rand_seqs[-1].append(''.join(np.random.choice(nucs,
                                           len(self.seqs[i]), p=probs)))

    def motif_significance(self, n=100, method='linear',
            fold_params={'win_size':None,'win_inc':None,
                         'sl_method':'loops'}):
        """Significance calculations wrapper.
        Args:
            n (int): Number of randomised sequences to calculate
            significance over.
            method (str): Method of significance calculation
            fold_params (dict): Parameters for folding 
        """

        # set parameters
        if not isinstance(n,int) or isinstance(n,bool):
            raise InputError('n must be integer!')
        self._nruns = n
        self._sig = True
        methods=['linear','ps']
         
        if not isinstance(method,str):
            raise InputError('method must be string')
        elif method not in methods:
            raise InputError('method not found')

        self._seq_randomiser()

        # run through calculations
        if method == 'linear':
            self._sig_calculator()
        if method == 'ps':
            self._fold_rand(fold_params)
            self.ps_sig_calculator(sl_method=fold_params['sl_method'])
        
        #output data
        self._sig_output(method)
    
    def _fold_rand(self,fold_params):
        """Fold the randomly generated sequences"""

        #iterate through the sequences
        for i in range(len(self.seqs)):

            #add reference sequence to the record list
            seq_recs = [SeqRecord(Seq(self.seqs[i]),
                        id = self.seq_names[i])]

            for j in range(len(self._rand_seqs[i])):

                 seq_recs.append(SeqRecord(Seq(self._rand_seqs[i][j]),
                      id = self.seq_names[i]+'_rand'+str(j)))
            self._fold(seq_recs,fold_params)

    def _fold(self,seq_recs,fold_params):
        """Folding algorithm
        Args:
            n (int): Number of randomised sequences to calculate
            significance over.
            seq_recs (:obj: 'list' of :obj: 'SeqRecord'): Record of sequences
            to be folded.
            fold_params (dict): Parameters for folding.
        """

        #create Twrap object for folding
        fold_pack = Twrap(seq_recs) 

        #fold sequences
        fold_pack.fold(win_size=fold_params['win_size'],
                       win_inc=fold_params['win_inc'],
                       sl_method=fold_params['sl_method'])
        
    def ps_sig_calculator(self,sl_method='loops'):
        """calculate significance of motif occuring in PS

        Args:
            sl_method (str): Stemloop determination method
        """

        #list methods for stemloop/loop analysis
        methods=['stemloops','loops']
        self._sig=True
        #make sure method selected is readable and correct
        if not isinstance(sl_method,str):
            raise InputError('sl method must be string!')
        elif sl_method not in methods:
            raise InputError('sl method not found!')
    
        #iterate through sequences
        for i in range(len(self.seqs)):

            #list stemloop files
            sl_files = glob.glob(sl_method+'_output/'+self.seq_names[i]+'*')
            #isolate reference file
            ref_file = [ii for ii in sl_files if ii == sl_method+\
                    '_output/' + self.seq_names[i] + '_'+sl_method+'.csv'][0]
            print(ref_file)
            #get reference value
            self._ps_mval.append(self._ps_motif_counter(ref_file))

            #initialise random values
            self._ps_mval_rand.append([])

            #get random seq files
            rand_files = [ii for ii in sl_files if ii != sl_method + \
                    '_output/' + self.seq_names[i] + '_'+sl_method+'.csv']

            #read through and calculate
            for j in rand_files:
                self._ps_mval_rand[i].append(self._ps_motif_counter(j))
                        
            #calculate random mean and standard dev
            self._ps_mean.append(np.mean(self._ps_mval_rand[i]))
            self._ps_std.append(np.std(self._ps_mval_rand[i]))

            # calculate Z-score of match count on investigated sequence
            if self._ps_std[i] != 0:

                self._ps_zscore.append(
                    (self._ps_mval[i] - self._ps_mean[i]) / self._ps_std[i])
            else:
                self._ps_zscore.append(0)

        # output significance values
        self._sig_output('ps')

    def _ps_motif_counter(self,file):
        """Read stemloop data files and extract motif occurences
        
        Args:
            file (str): File containing stemloop data
        """

        #open datafile
        data = pd.read_csv(file)
        data2=self._sl_data_reader(data)
        #initialise counter
        counter = 0
        freq_col = data.columns.get_loc('freq')
        #return counter
        return np.sum(data2[:,int(freq_col)])

    def _ps_motif_filter(self):
        """Read stemloop data files and extract motif occurences
        
        Args:
            file (str): File containing stemloop data
        """
        self._sig = True
        #open datafile
        files_in = glob.glob('*loops_output/*')

        files_in = [ff for ff in files_in if 'loops.csv' in ff]
        for file in files_in:
            data = pd.read_csv(file)

            data2 = self._sl_data_reader(data)
            np.savetxt(file.split('.csv')[0] + '_' + self._motif+'.csv',
            data2.astype(str),delimiter=',',
            fmt = '%s',header=','.join(data.columns))

    def _sl_data_reader(self,data):
        """Read stemloop data files and extract motif occurences
        
        Args:
            file (str): File containing stemloop data
        """
        i = 0
        #initialise counter
        #read through loop data
        data2 = []
        for i in range(len(data)):
            self._temp_match=[]
            #scan each loop sequence for motif occurence
            self._mscan(data['loop_seq'][i])
            #if motif in loop count frequency
            if  len(self._temp_match)>0:
                if len(data2) == 0:
                    data2 = data.values[i,:]
                else:
                    data2 = np.vstack((data2,data.values[i,:]))

        return data2           
 
    def _sig_calculator(self):
        """Calculate the significance of motif occurence
        """

        # iterate through sequences
        for i in range(len(self.seqs)):

            # set actual value
            if self.match[i].ndim == 1:
                self._aval.append(1)
            else:
                self._aval.append(len(self.match[i]))

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

    def _sig_output(self,method='linear'):
        """Writes significance values to file

        Args:
            method (str): Significance calculation method
        """
        methods=['linear','ps']
        #check method input correct
        if not isinstance(method,str):
            raise InputError('method must be string')
        elif method not in methods:
            raise InputError('method not found')

        #based on method set file name
        if method == 'ps':
            file_name=self._motif + '_ps_sig.csv'
            
            file_out = open(file_name, 'w')

            # write headers
            file_out.write('seq,count,mean,std,zscore\n')

            # iterate through significance values and write to file
            for i in range(len(self.seqs)):
                file_out.write(
                    '%s,%d,%f,%f,%f\n' %
                    (self.seq_names[i],
                     self._ps_mval[i],
                        self._ps_mean[i],
                        self._ps_std[i],
                        self._ps_zscore[i]))

            file_out.close()
        else:
            file_name=self._motif + '_sig.csv'

            file_out = open(file_name, 'w')

            # write headers
            file_out.write('seq,count,mean,std,zscore\n')

            # iterate through significance values and write to file
            for i in range(len(self.seqs)):
                file_out.write(
                    '%s,%d,%f,%f,%f\n' %
                    (self.seq_names[i],
                     self._aval[i],
                        self._mean[i],
                        self._std[i],
                        self._zscore[i]))

            file_out.close()

    def _data_output(self):
        """Output motif matches to file.
        """

        # iterate through sequences
        for j in range(len(self.seq_names)):

            # open file for each sequence
            file_out = open(
                self.seq_names[j] + '_' + self._motif + '_occurences.csv', 'w')

            # set header
            file_out.write('position,motif\n')

            # write matches to file

            if self.match[j].ndim == 1:
                file_out.write('%d,%s\n' %
                               (int(self.match[j][0]), self.match[j][1]))
            else:
                for i in range(len(self.match[j])):
                    file_out.write('%d,%s\n' %
                                   (int(self.match[j][i, 0]), self.match[j][i, 1]))

            file_out.close()

    def stats_output(self):
        print('hello')
        for j in range(len(self.seq_names)):
            file_out = open(self.seq_names[j] + '_statistics.csv', 'w')
            file_out.write('metric,value\n')
            file_out.write('length,%d\n' % (len(self.seqs[j])))
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

            file_out.close()
