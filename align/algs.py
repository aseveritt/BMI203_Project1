
import pandas as pd
import numpy as np
import sys


##### -------------------------------------------------------------------------------  #####

#Functions involving fasta

class FastaRecord(object):
	'''Stores name and sequence of single fasta record
    
    Parameters:
        sequence file (str): FASTA file located in sequence directory to be loaded

    Attributes:
        seq (str): capitalized protein sequence which only contans valid characters
        name (str): read or protein name
	'''
	def __init__(self, sequencefile):
		self.title, self.seq = read_fasta(sequencefile)

def checkSeq(seq):
	'''Helper function thats flexible to check protein sequences are acceptible. 
	
	| Currently checks if sequence contains only alphabetical letters, 
	| formats into upper case only, 
	| and checks if only protein AA bases are given. (with x ans z additonally added)

	Parameters: 
		seq (str) : protein sequence
	
	Returns: 
		uppercase of seqeunce provided or error
    '''

	ok = dict.fromkeys("ARNDCQEGHILKMFPSTWYVXZ") #adding X and Z but not happy about it. 

	if seq.isalpha():
		seq = seq.upper()
		if all(c in ok for c in seq):
			return(seq)
		else: print("You've entered an invalid AA base. Please try again."); sys.exit(1)
	else: print("You've entered a sequence with non-alphabetical letters. Please try again."); sys.exit(1)


def read_fasta(sequencefile):
	'''Reads in a single fasta sequnce from a file which must be in the sequences directory. 
	
	| Not adapted for other directories or having multiple sequences in one file. 
	| Not adapted for .fastq

	Parameters: 
		sequencefile (str) : fasta file to read in. 
	
	:returns:
		- name (str): read name
		- sequence (str) : upper case, AA acid sequence
    '''
	with open(sequencefile) as fp: 
		name, seq = None, []
		for line in fp:
			line = line.rstrip()
			if line.startswith(">"):
				name, seq = line, []
			else:
				seq.append(line)
	if name: 
		final_seq = ''.join(seq)
		return name, checkSeq(final_seq)
	else : print("Check File is not empty"); sys.exit(1)

##### -------------------------------------------------------------------------------  #####


#Helpers for unit testing

def generate_DFcmd_for_pytest(df, pandas_as='pd'):
    """
    Usage: 
    a = PairwiseAligner(method="SW", gap_cost=5, substitutionMatrix = "BLOSUM62")
	df = a.sub_mat
	print(gencmd(df))

	Takes an input df and outputs the command neccesssary to recreate it so I can use it for unit testing. 
    """

    #Credit here : https://stackoverflow.com/questions/24005466/given-a-pandas-dataframe-is-there-an-easy-way-to-print-out-a-command-to-generat
    from types import MethodType
    from pandas import DataFrame, MultiIndex

    if pandas_as:
    	pandas_as += '.'
    index_cmd = df.index.__class__.__name__
    if type(df.index)==MultiIndex:
    	index_cmd += '.from_tuples({0}, names={1})'.format([i for i in df.index], df.index.names)
    else:
    	index_cmd += "({0}, name='{1}')".format([i for i in df.index], df.index.name)
    	return 'DataFrame({0}, index={1}{2}, columns={3})'.format([[xx for xx in x] for x in df.values],
                                                                pandas_as,
                                                                index_cmd,
                                                                [c for c in df.columns])




##### -------------------------------------------------------------------------------  #####



def load_substitution_matrix(subchoice):
	'''Loads in substitution matrix of choice. 

	| From the four possible options that are located in the scoring matrices directory,
	| read them in as a matrix. I choose to do this way so its more robust towards errors in syntax. 
	| Then convert to a pandas dataframe to make indexing easier. 

	Parameters: 
		subchoice (str): passed from PairwiseAligner class. must be either: PAM250, PAM100, BLOSUM50, BLOSUM62

	Returns: 
		scoring matrix (df) : dataframe where rows and column are AA index
    '''

	if subchoice in set(["PAM250","PAM100", "BLOSUM50", "BLOSUM62"]):
		myfile = "scoring_matrices/" + subchoice + ".mat"
		with open(myfile, 'r') as f:
			score_mat = [[num for num in line.strip().split(' ') if num!= ""] for line in f if not "#" in line]
		df = pd.DataFrame(score_mat)
		new_indexes = df.iloc[0]; df = df[1:]
		df.columns = new_indexes; df.index = new_indexes
		return(df)
	else: print("Please select a substitution matrix from the following list: PAM250, PAM100, BLOSUM50, BLOSUM62"); sys.exit(1)



def print_like_blast(aligned_seq1, aligned_seq2, print_flag):
	'''
	Create an output format similar to blast. 
	Print percent identities as well as percent gaps. 

	Parameters: 
		aligned_seq1 (str): string representing the aligned segments with gaps where neccessary. passed from backtrace method. 
		aligned_seq2 (str): string representing the aligned segments with gaps where neccessary. passed from backtrace method. 
		print_flag (bool): whether to print output to screen (True) or just return the number of identities and gaps (False)

	:returns:
		- idents (int): raw number of positive identities
		- gaps (int): number of gaps in aligned sequences
	'''
	#CREDIT: https://gist.github.com/radaniba/11019717 


	idents, gaps, mismatches = 0, 0, 0
	alignment_string = []
	for base1, base2 in zip(aligned_seq1, aligned_seq2):
		if base1 == base2:
			alignment_string.append('|')
			idents += 1
		elif '-' in (base1, base2):
			alignment_string.append(' ')
			gaps += 1
		else:
			alignment_string.append(':')
			mismatches += 1
	alength = len(aligned_seq1)
	alignment_str = ''.join(alignment_string)

	if not print_flag:
		return idents, gaps

	print()
	print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents, alength, idents / alength, gaps, alength, gaps / alength))
	print()
	for i in range(0, alength, 60):
		seq1_slice = aligned_seq1[i:i+60]
		print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
		print('             {0}'.format(alignment_str[i:i+60]))
		seq2_slice = aligned_seq2[i:i+60]
		print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
		print()
	return idents, gaps

##### -------------------------------------------------------------------------------  #####


##### -------------------------------------------------------------------------------  #####


class PairwiseAligner:
	"""Performs a Pairwise Alignment with of two provided sequences.
	
	Currently adapted for Smith-Waterman with affine gap penalties and Needleman Wunsch with affine gap penalties. 

	| Obtains helper functions:
	| initialize_matrices(): nothing, passes to child classes
	| init_v(): initialize the pointer matrix. same between two classes. 
	| match(): score of match/mismatch
	| fill_matrices(): fill all 3 matrices in nested for loop. 
	| backtrace(): trace through pointer matrix. 
	| number_matches():
	| number_gaps():

     Attributes:
        gap_open (float): gap opening cost. Suggested values fall between -1 and -20
        gap_extension (float): gap extension cost. Suggested values fall between -1 and -5
        substitutionMatrix (str): choice of four substiution matrices "BLOSUM62" (default), "BLOSUM50" ,"PAM250", "PAM100"
    """

	def __init__(self, gap_open=1, gap_extension = 1, substitutionMatrix = "BLOSUM62"):
		self.gap_open = gap_open
		self.gap_ext = gap_extension
		self.sub_mat = load_substitution_matrix(substitutionMatrix)

	def initialize_matrices(self):
		#just pass to child
		pass

	def init_v(self, i, j):
		''' Initialize Pointer Matrix for traceback. 
		
		| v[0,0] set to "end"
		| top row set to "left"
		| leftmost column set to "up"
		| all others set to "tbd" which will be overwritten by distance matrix function. 

		| No difference between global and local. Internal function called by distance matrix. 

		Parameters: 
			i, j (int,int): position in matrix to fill. 
		'''
		if j == 0 and i == 0:  return "end"
		elif i == 0 and j > 0: return "left"
		elif j == 0 and i > 0: return "up"
		else: return "tbd"

	def match(self, charA, charB, method="real", match=1, mismatch=-1):
		''' Match score between two AAs
	
		| Using a substitution matrix (method="real"), return the score of pairing the two AA together. 
		| Using a priori values (method = "testing"), return either the match or mismatch score. 

		Parameters: 
			charA (str): AA character (must be uppercase and found in sub_mat)
			charB (str): AA character (must be uppercase and found in sub_mat)
			method (str): either "real" (default) or "testing" depending on use of substitution matrix
			match (int): default 1 used if method == "testing"
			mismatch (int): default -1 used if method == "testing"

		Returns:
			(int)
		'''
		if method=="real":
			return int(self.sub_mat.loc[charA, charB])
		if method == "testing" :
			if charA == charB: return  match
			else: return mismatch

	def fill_matrices(self, method="real"):
		''' Builld all scoring and traceback matrices needed for performing a pairwise alignment. 
	
		| Four matrices are initialized of size len(seqB) + 1 * len(seqA) + 1. 

		| For all combination of sequence AAs, the value of a possible  match/mismatch, gap in seqA, or gap in seqB is obtained. 
		| In NW, the max value of these three possible options is saved in the total M matrix.
		| In SW, the max value of these three possible options and zero  (as to not go negtive) is saved in the total M matrix.    
		| Depending on which matrix the top score originates from, the 'val' matrix is filled to determine what direction a user will go in the backtrace. 
	
		| As the M matrix is filled, the value and location of the cell with the highest value is stored. 

		Parameters: 
			self (PairwiseAligner): requires variables gap_open, gap_ext, seqA, seqB, and method be set. 
			method (str): either "real" (default) or "testing" depending on if substitution matrix is being used. 

		:returns:
			- val (matrix): traceback matrix containing cells of either ("diag", "up","left","end") as instructions for backtrace
			- M (matrix): Matrix representing the maxiumum cell obtained by either a (gap in either orientation X[i,j] or Y[i,j] or a match/mismatch score at i,j
			- max_val (int): the maximum value observed in the M matrix -- this is the raw alignment score in most situations. 
			- max_pos (str): row and column in M matrix where the maximum value is observed. takes the form (5,6)
		'''

		#I know this is technically not using code efficiently, but so much easier to read.
		gap_open = self.gap_open
		gap_ext = self.gap_ext
		M = self.M
		X = self.X
		Y = self.Y
		val = self.TB

		max_val = 0
		max_pos = None
		#print(M); print(X); print(Y); print(val)

		for i in range(1, len(self.seqB) + 1):
			for j in range(1, len(self.seqA) + 1):
				tmp = self.match(self.seqA[j - 1], self.seqB[i - 1], method) + M[i-1][j-1] #tmp is match/mismatch
				Y[i][j] = max(
						  (M[i][j-1] + gap_open), 
			 			  (Y[i][j-1] + gap_ext))
				X[i][j] = max(
						  (M[i-1][j] + gap_open), 
						  (X[i-1][j] + gap_ext)
						  )

				if self.flag == "NW" : ugh = [tmp, X[i][j], Y[i][j]]
				if self.flag == "SW" : ugh = [tmp, X[i][j], Y[i][j], 0]

				M[i][j]  = max(ugh)	
				if M[i][j] > max_val:  				  #calculate it so that we know where to start traceback
					max_val = M[i][j]
					max_pos = (i,j)

				maxmat = ugh.index(max(ugh))			#first occurence given if tie, so ties will go to diagnols.
				if maxmat == 0:   val[i][j] = "diag"
				elif maxmat == 1: val[i][j] = "up"
				elif maxmat == 2: val[i][j] = "left"
				elif maxmat == 3: val[i][j] = "end"				#will never reach in NW
				else: print("ugh bro howd this happen"); exit()

		if method == "testing" : return [M]
		return [val, M, max_val, max_pos]

	def backtrace(self, val):
		''' Traceback through alignment matrix to find structure of top pairwise of provided sequences.

		Parameters: 
			val (matrix): traceback matrix filled with either "diag","up","left","end" passed from  distance_matrix function

		Returns:
			aligned_sequences (str): string representing the optimal alignment of the provided sequences and traceback matrix. eg. ('AAT-','AATC')
		'''

		i,j = self.top_pos
		sequ1 = ''
		sequ2 = '' 

		while (i>0 or j>0):
			currval = val[i][j]
			if currval == "diag": #diag
				sequ1 += self.seqA[j-1]
				sequ2 += self.seqB[i-1]
				i -= 1; j -= 1
			elif currval == "up": #down
				sequ1 += '-'
				sequ2 += self.seqB[i-1]
				i -= 1
			elif currval == "left": #down
				sequ1 += self.seqA[j-1]
				sequ2 += '-'
				j -= 1
			elif currval == "end":
				break
			else:
				print("How tf did we get here")
				print(currval); print(i, j)
				break

		if len(sequ1) != len(sequ2):
			raise ValueError('issue in traceback: alignments not the same length')

		aligned_sequences = ''.join(reversed(sequ1)), ''.join(reversed(sequ2))
		return aligned_sequences

	def align(self, seqA, seqB, print_flag=True):
		''' Align two sequences together

		Parameters: 
			self (PairwiseAligner): with attributes gap_open, gap_ext, sub_mat, and method defined.
			seqA (FastaRecord): Fasta Record representing sequence of AA characters and name
			seqB (FastaRecord): Fasta Record representing sequence of AA characters and name
			print_flag (bool): default (True) will print alignment to screen

		Returns:
			self (PairwiseAligner): self with new attributes num_matches, num_gaps, top_pos, and top_score added
		'''

		self.seqA = seqA.seq
		self.seqB = seqB.seq
		self.initialize_matrices()
		traceback_mat, M_mat, score, position = self.fill_matrices()

		if self.flag == "SW": 
			self.top_pos = position
			self.top_score = score

		if self.flag == "NW": 
			self.top_pos = (len(self.seqB), len(self.seqA))
			self.top_score = M_mat[len(self.seqB)][len(self.seqA)]

		seq1_aligned, seq2_aligned = self.backtrace(traceback_mat)
		match, gaps = print_like_blast(seq1_aligned, seq2_aligned, print_flag)
		self.num_matches = match
		self.num_gaps = gaps
		return(self)

	def raw_score(self):
		return self.top_score

	def number_matches(self):
		return self.num_matches

	def number_gaps(self):
		return self.num_gaps


class SmithWaterman(PairwiseAligner):
	"""Computes Smith-Waterman measure.

    | The Smith-Waterman algorithm with affine gap pentalty performs local sequence alignment between two strings. 
    | This will not align the entire sequence length, but the site with highest similarity. 

    Parameters:
        gap_open (float): gap opening cost. Suggested values fall between -1 and -20
        gap_extension (float): gap extension cost. Suggested values fall between -1 and -5
        substitutionMatrix (str): choice of four substiution matrices "BLOSUM62" (default), "BLOSUM50" ,"PAM250", "PAM100"
	"""

	def __init__(self, gap_open, gap_extension, substitutionMatrix, method="SW"):
		self.flag = method
		PairwiseAligner.__init__(self, gap_open, gap_extension, substitutionMatrix)
	
	def init_m(self, i, j):
		if j == 0 and i == 0:  return 0  #M[0,0]
		if i == 0 and j > 0: return 0 
		if j == 0 and i > 0: return 0
		else: return -float("inf")

	def init_y(self, i, j):
		if i == 0 and j == 0: return self.gap_open
		else: return -float("inf")

	def init_x(self, i, j):
		if i == 0 and j == 0: return self.gap_open
		else: return -float("inf")

	def initialize_matrices(self):
		dim_i = len(self.seqB) + 1
		dim_j = len(self.seqA) + 1

		self.X = [[self.init_x(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		self.Y = [[self.init_y(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		self.M = [[self.init_m(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		self.TB = [[self.init_v(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		return (self)


class NeedlemanWunsch(PairwiseAligner):
	"""Computes Needleman-Wunsch measure.

    The Needleman-Wunsch algorithm with affine gap pentalty performs global sequence alignment between two strings. 
    This will align the entire sequence length.  

    Parameters:
        gap_open (float): gap opening cost. Suggested values fall between -1 and -20
        gap_extension (float): gap extension cost. Suggested values fall between -1 and -5
        substitutionMatrix (str): choice of four substiution matrices "BLOSUM62" (default), "BLOSUM50" ,"PAM250", "PAM100"
	"""

	def __init__(self, gap_open, gap_extension, substitutionMatrix, method="NW"):
		self.flag = method
		PairwiseAligner.__init__(self, gap_open, gap_extension, substitutionMatrix)
	
	def init_m(self, i, j):
		if j == 0 and i == 0:  return 0  #M[0,0]
		else: return -float("inf")

	def init_y(self, i, j):
		if i == 0 and j == 0:  return self.gap_open
		elif i == 0 and j > 0: return (self.gap_open + (self.gap_ext * j))
		else: return -float("inf")
	def init_x(self, i, j):
		if i == 0 and j == 0:  return self.gap_open
		elif j == 0 and i > 0: return (self.gap_open + (self.gap_ext * i))
		else: return -float("inf")

	def initialize_matrices(self):
		dim_i = len(self.seqB) + 1
		dim_j = len(self.seqA) + 1

		self.X = [[self.init_x(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		self.Y = [[self.init_y(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		self.M = [[self.init_m(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		self.TB = [[self.init_v(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
		return (self)




##### -------------------------------------------------------------------------------  #####



if __name__ == "__main__":
	nw = NeedlemanWunsch(gap_open=-10, gap_extension = -0.5, substitutionMatrix = "BLOSUM50")
	sw = SmithWaterman(gap_open=-10, gap_extension = -1, substitutionMatrix = "BLOSUM62")

	b= sw.align(FastaRecord("sequences/prot-0004.fa"), FastaRecord("sequences/prot-0004.fa"), print_flag=True)
	print(b.raw_score())
	print(b.number_matches())
	print(b.number_gaps())
	b= nw.align(FastaRecord("sequences/prot-0004.fa"), FastaRecord("sequences/prot-0004.fa"), print_flag=True)
	print(b.raw_score())
	print(b.number_matches())
	print(b.number_gaps())

	nw = NeedlemanWunsch(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM62")
	b = sw.align(FastaRecord("sequences/test5.fa"), FastaRecord("sequences/test6.fa"), print_flag=True)
	print(b.top_pos)



	
