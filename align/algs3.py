
import pandas as pd
import numpy as np
import sys


##### -------------------------------------------------------------------------------  #####

#Functions involving fasta

class FastaRecord(object):
    def __init__(self, sequencefile):
        self.title, self.seq = read_fasta(sequencefile)

def checkSeq(seq):
	'''Helper function thats flexible to check protein sequences are acceptible. 
	Currently checks if sequence contains only alphabetical letters, 
	formats into upper case only, 
	and checks if only protein AA bases are given. 

	Args: seq (str) : protein sequence
	Output: uppercase of seqeunce provided or error
    '''

	ok = dict.fromkeys("ARNDCQEGHILKMFPSTWYV")

	if seq.isalpha():
		seq = seq.upper()
		if all(c in ok for c in seq):
			return(seq)
		else: print("You've entered an invalid AA base. Please try again."); sys.exit(1)
	else: print("You've entered a sequence with non-alphabetical letters. Please try again."); sys.exit(1)


def read_fasta(sequencefile):
	'''Reads in a single fasta sequnces from a file which must be in the sequences directory. 
	Not adapted for other directories or having multiple sequences in one file. 
	Not adapted for .fastq

	Args: 
		sequencefile (str) : fasta file to read in. 
	Output: 
		name (str): read name
		sequence (str) : upper case, AA acid sequence
    '''
	#with open('../' + sequencefile) as fp: 
	with open(sequencefile) as fp: 
		name, seq = None, []
		for line in fp:
			line = line.rstrip()
			if line.startswith(">"):
				#if name: yield (name, ''.join(seq))
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
	Credit here : https://stackoverflow.com/questions/24005466/given-a-pandas-dataframe-is-there-an-easy-way-to-print-out-a-command-to-generat
    """
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

	From the four possible options that are located in the scoring matrices directory,
	read them in as a matrix. I choose to do this way so its more robust towards errors in syntax. 
	Then convert to a pandas dataframe to make indexing easier. 

	Args: subchoice (str): passed from PairwiseAligner class. 

	Attributes: scoreing matrix (df) : dataframe where rows and column are AA index
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




def build_scoring_matrix(self, seqA, seqB):

	#initialize scoring matrix
	#SW first row and column at zero.
	scoring_mat = np.zeros((len(seqA) + 1, len(seqB) + 1), dtype=np.float) #same in both
	scoring_mat[0][1] = - self.gap_open
	scoring_mat[1][0] = - self.gap_open

	print(scoring_mat)
	exit()
	max_value = 0

	if self.method=="NW" : max_pos = (len(seqA) , len(seqB) ) #global means traceback always starts at bottom corner
	else: max_pos = None

	for i in range(1, len(seqA) + 1): # 1 indexed
		for j in range(1, len(seqB) + 1):
			match = scoring_mat[i - 1, j - 1] + int(self.sub_mat.loc[seqA[i - 1], seqB[j - 1]]) 
			delete = scoring_mat[i - 1, j] - self.gap_open
			insert = scoring_mat[i, j - 1] - self.gap_open
			if self.method=="NW" : 
				scoring_mat[i, j] = max(match, delete, insert)
				max_value = max(max_value, scoring_mat[i, j])
			elif self.method=="SW": 
				scoring_mat[i, j] = max(0, match, delete, insert) #add zero so it doesn't go negative. 
				if scoring_mat[i, j] > max_value:  				  #calculate it so that we know where to start traceback
					max_value = scoring_mat[i, j]
					max_pos = (i,j)
			else:
				print("ERROR in scoring matrix"); exit()

	return scoring_mat, max_value, max_pos



def traceback(flag, score_matrix, start_pos, seq1, seq2):

    '''Find the optimal path through the matrix.
    This function traces a path from the bottom-right to the top-left corner of
    the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
    or both of the sequences being aligned. Moves are determined by the score of
    three adjacent squares: the upper square, the left square, and the diagonal
    upper-left square.
    WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 1
        left:     gap in sequence 2
    '''

    aligned_seq1 = []
    aligned_seq2 = []
    x, y         = start_pos
    move         = next_move(flag, score_matrix, x, y)
    #print((x,y))
    #print(move)
    while move != "END":
        if move == "DIAG":
        	#if diagnol, subtract 1 from both coordinates and add to aligned seq. 
            aligned_seq1.append(seq1[x - 1]) 
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == "UP":
        	# subtract from columns, but not row. 
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
        else:
        	#Move left, so subtract from rows
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1
        #print((x,y))
        move = next_move(flag, score_matrix, x, y) #pick next move until we reach the end. 
        #print(move)


    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq2[y - 1])
    a = ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))
    
    if len(aligned_seq1) != len(aligned_seq2):
    	raise ValueError('issue in traceback: alignments not the same length')

    return a

def next_move(flag, score_matrix, x, y):
	if x == 1 and y == 1: return "END"

	diag = score_matrix[x - 1][y - 1]
	up   = score_matrix[x - 1][y]
	left = score_matrix[x][y - 1]
	if diag >= up and diag >= left:     # Tie goes to the DIAG move.
		if flag == "SW": return "DIAG" if diag != 0 else "END"  
		elif flag == "NW" : return "DIAG"
	elif up > diag and up >= left:      
		#return "UP" if up != 0 else "END"  
		if flag == "SW": return "UP" if up != 0 else "END"  
		elif flag == "NW" : return "UP"    
	elif left > diag and left > up:
		#return "LEFT" if left != 0 else "END"   
		if flag == "SW": return "LEFT" if left != 0 else "END"  
		elif flag == "NW" : return "LEFT"
	else:
        # Execution should not reach here.
		print("WHY ARE WE HERE")
		exit()


def print_like_blast(aligned_seq1, aligned_seq2):
	'''
	STOLEN: AMANDA CREDIT
	https://gist.github.com/radaniba/11019717
	Construct a special string showing identities, gaps, and mismatches.
	This string is printed between the two aligned sequences and shows the
	identities (|), gaps (-), and mismatches (:). As the string is constructed,
	it also counts number of identities, gaps, and mismatches and returns the
	counts along with the alignment string.
	AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
	::||::::::||:|::::::: |:  :||:|   <-- alignment string
	CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
	'''
	# Build the string as a list of characters to avoid costly string
	# concatenations

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
	return 


##### -------------------------------------------------------------------------------  #####


class PairwiseAligner:
	def __init__(self, method= "SW", gap_open=1, gap_extension = 1, substitutionMatrix = "PAM250"):
		self.gap_open = gap_open
		self.gap_ext = gap_extension
		self.method = method
		self.sub_mat = load_substitution_matrix(substitutionMatrix)
	
	def align(self, seqA, seqB, print_flag):
		#check inputs
		seqA = seqA.seq
		seqB = seqB.seq
		score_matrix, max_score, max_pos = build_scoring_matrix(self, seqA, seqB)
		#print(max_pos)
		#print(score_matrix)
		self.highest_score = max_score
		seq1_aligned, seq2_aligned = traceback(self.method, score_matrix, max_pos, seqA, seqB)
		#print(seq1_aligned)
		if print_flag: print_like_blast(seq1_aligned, seq2_aligned)
		return(self)

#	def alignNW(self, seqA, seqB):
#		seqA = checkSeq(seqA)
#		seqB = checkSeq(seqB)
#		score_matrix, max_score, max_pos = build_scoring_matrix(self, seqA, seqB)
#		self.highest_score = max_score
#		max_pos = (len(seqA) + 1, len(seqB) + 1)
#		seq1_aligned, seq2_aligned = traceback(self.method, score_matrix, max_pos, seqA, seqB)
#		print_like_blast(seq1_aligned, seq2_aligned)
#		return(self)

	def get_top_score():
		return self.highest_score

#	def get_top_match():
#		pass

class SmithWaterman(PairwiseAligner):

	"""Computes Smith-Waterman measure.

    The Smith-Waterman algorithm performs local sequence alignment; that is, for determining similar regions
    between two strings. Instead of looking at the total sequence, the Smith–Waterman algorithm compares segments of
    all possible lengths and optimizes the similarity measure. See the string matching chapter in the DI book (Principles of Data Integration). 

    Args:
        gap_cost (float): Cost of gap (defaults to 1.0).
        sim_func (function): Similarity function to give a score for the correspondence between the characters (defaults
                             to an identity function, which returns 1 if the two characters are the same and 0 otherwise).

    Attributes:
        gap_cost (float): An attribute to store the gap cost.
        sim_func (function): An attribute to store the similarity function.
    """
    
	#https://anhaidgroup.github.io/py_stringmatching/v0.3.x/_modules/py_stringmatching/similarity_measure/smith_waterman.html#SmithWaterman
	#https://anhaidgroup.github.io/py_stringmatching/v0.3.x/_modules/py_stringmatching/similarity_measure/affine.html#Affine
	pass

class NeedlemanWunsch(PairwiseAligner):
	#https://anhaidgroup.github.io/py_stringmatching/v0.3.x/_modules/py_stringmatching/similarity_measure/needleman_wunsch.html#NeedlemanWunsch
	pass








if __name__ == "__main__":

	a = PairwiseAligner(method="SW", gap_open=5, gap_extension = 1, substitutionMatrix = "BLOSUM62")
	#df = a.sub_mat
	print(a.align(FastaRecord("sequences/prot-0004.fa"), FastaRecord("sequences/prot-0008.fa"), print_flag=True).highest_score)

	a = PairwiseAligner(method="NW", gap_open=5, substitutionMatrix = "BLOSUM62")
	#df = a.sub_mat
	print(a.align(FastaRecord("sequences/prot-0004.fa"), FastaRecord("sequences/prot-0008.fa"), print_flag=True).highest_score)

