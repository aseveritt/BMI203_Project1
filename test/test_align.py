import pytest
from align import algs
import pandas as pd

@pytest.fixture
def some_relevant_data():
	return np.ones(10)

def test_fasta_io():
	#what check sequnces is supposed to catch
	assert algs.checkSeq("AAAATTARRYY") == "AAAATTARRYY", "Failing base case" #fixes nothing
	assert algs.checkSeq("aaaATTARRYY") == "AAAATTARRYY", "Failing to captialize" 
	with pytest.raises(SystemExit):
		algs.checkSeq("aaa ATTARRYY") #should fail if non alpha numeric
	with pytest.raises(SystemExit):
		algs.checkSeq("UaaaATTARRYY") #should fail if non-AA bases appear


	#what read_fasta is supposed to catch
	assert algs.read_fasta("sequences/prot-0047.fa")  == ('>d2int__ 1.22.1.2.2', 'HKCDITLQEIIKTLNSLTEQKTLCTELTVTDIFAASKNTTEKETFCRAATVLRQFYSHHEKDTRCLGATAQQFHRHKQLIRFLKRLDRNLWGLAGLNSCPVKEANQSTLENFLERLKTIMREKYSKCSS'), "Failing base case"
	with pytest.raises(FileNotFoundError):
		algs.read_fasta("sequences/imnothere.fa")

def test_sub_matrix_io():
	correct_df = pd.DataFrame([['4', '-1', '-2', '-2', '0', '-1', '-1', '0', '-2', '-1', '-1', '-1', '-1', '-2', '-1', '1', '0', '-3', '-2', '0', '-2', '-1', '0', '-4'], 
		['-1', '5', '0', '-2', '-3', '1', '0', '-2', '0', '-3', '-2', '2', '-1', '-3', '-2', '-1', '-1', '-3', '-2', '-3', '-1', '0', '-1', '-4'], 
		['-2', '0', '6', '1', '-3', '0', '0', '0', '1', '-3', '-3', '0', '-2', '-3', '-2', '1', '0', '-4', '-2', '-3', '3', '0', '-1', '-4'], 
		['-2', '-2', '1', '6', '-3', '0', '2', '-1', '-1', '-3', '-4', '-1', '-3', '-3', '-1', '0', '-1', '-4', '-3', '-3', '4', '1', '-1', '-4'], ['0', '-3', '-3', '-3', '9', '-3', '-4', '-3', '-3', '-1', '-1', '-3', '-1', '-2', '-3', '-1', '-1', '-2', '-2', '-1', '-3', '-3', '-2', '-4'], 
		['-1', '1', '0', '0', '-3', '5', '2', '-2', '0', '-3', '-2', '1', '0', '-3', '-1', '0', '-1', '-2', '-1', '-2', '0', '3', '-1', '-4'], ['-1', '0', '0', '2', '-4', '2', '5', '-2', '0', '-3', '-3', '1', '-2', '-3', '-1', '0', '-1', '-3', '-2', '-2', '1', '4', '-1', '-4'], 
		['0', '-2', '0', '-1', '-3', '-2', '-2', '6', '-2', '-4', '-4', '-2', '-3', '-3', '-2', '0', '-2', '-2', '-3', '-3', '-1', '-2', '-1', '-4'], ['-2', '0', '1', '-1', '-3', '0', '0', '-2', '8', '-3', '-3', '-1', '-2', '-1', '-2', '-1', '-2', '-2', '2', '-3', '0', '0', '-1', '-4'], 
		['-1', '-3', '-3', '-3', '-1', '-3', '-3', '-4', '-3', '4', '2', '-3', '1', '0', '-3', '-2', '-1', '-3', '-1', '3', '-3', '-3', '-1', '-4'], ['-1', '-2', '-3', '-4', '-1', '-2', '-3', '-4', '-3', '2', '4', '-2', '2', '0', '-3', '-2', '-1', '-2', '-1', '1', '-4', '-3', '-1', '-4'], 
		['-1', '2', '0', '-1', '-3', '1', '1', '-2', '-1', '-3', '-2', '5', '-1', '-3', '-1', '0', '-1', '-3', '-2', '-2', '0', '1', '-1', '-4'], ['-1', '-1', '-2', '-3', '-1', '0', '-2', '-3', '-2', '1', '2', '-1', '5', '0', '-2', '-1', '-1', '-1', '-1', '1', '-3', '-1', '-1', '-4'], 
		['-2', '-3', '-3', '-3', '-2', '-3', '-3', '-3', '-1', '0', '0', '-3', '0', '6', '-4', '-2', '-2', '1', '3', '-1', '-3', '-3', '-1', '-4'], ['-1', '-2', '-2', '-1', '-3', '-1', '-1', '-2', '-2', '-3', '-3', '-1', '-2', '-4', '7', '-1', '-1', '-4', '-3', '-2', '-2', '-1', '-2', '-4'], 
		['1', '-1', '1', '0', '-1', '0', '0', '0', '-1', '-2', '-2', '0', '-1', '-2', '-1', '4', '1', '-3', '-2', '-2', '0', '0', '0', '-4'], ['0', '-1', '0', '-1', '-1', '-1', '-1', '-2', '-2', '-1', '-1', '-1', '-1', '-2', '-1', '1', '5', '-2', '-2', '0', '-1', '-1', '0', '-4'], 
		['-3', '-3', '-4', '-4', '-2', '-2', '-3', '-2', '-2', '-3', '-2', '-3', '-1', '1', '-4', '-3', '-2', '11', '2', '-3', '-4', '-3', '-2', '-4'], ['-2', '-2', '-2', '-3', '-2', '-1', '-2', '-3', '2', '-1', '-1', '-2', '-1', '3', '-3', '-2', '-2', '2', '7', '-1', '-3', '-2', '-1', '-4'], 
		['0', '-3', '-3', '-3', '-1', '-2', '-2', '-3', '-3', '3', '1', '-2', '1', '-1', '-2', '-2', '0', '-3', '-1', '4', '-3', '-2', '-1', '-4'], ['-2', '-1', '3', '4', '-3', '0', '1', '-1', '0', '-3', '-4', '0', '-3', '-3', '-2', '0', '-1', '-4', '-3', '-3', '4', '1', '-1', '-4'], 
		['-1', '0', '0', '1', '-3', '3', '4', '-2', '0', '-3', '-3', '1', '-1', '-3', '-1', '0', '-1', '-3', '-2', '-2', '1', '4', '-1', '-4'], ['0', '-1', '-1', '-1', '-2', '-1', '-1', '-1', '-1', '-1', '-1', '-1', '-1', '-1', '-2', '0', '0', '-2', '-1', '-1', '-1', '-1', '-1', '-4'], 
		['-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '-4', '1']], index=pd.Index(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'], name='0'), 
		columns=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'])
	
	assert correct_df.equals(algs.load_substitution_matrix("BLOSUM62")) == True, "Failing base case" #check one of our knowns
	
	with pytest.raises(SystemExit):
		algs.load_substitution_matrix("doesntexist") #check if fails properly

def test_scoring_matrix_io():
	#initialization
	MIN=-float("inf")
	correct_X = [[-3, MIN, MIN, MIN, MIN, MIN], [-4, MIN, MIN, MIN, MIN, MIN], [-5, MIN, MIN, MIN, MIN, MIN], [-6, MIN, MIN, MIN, MIN, MIN]]
	correct_Y = [[-3, -4, -5, -6, -7, -8], [MIN, MIN, MIN, MIN, MIN, MIN], [MIN, MIN, MIN, MIN, MIN, MIN], [MIN, MIN, MIN, MIN, MIN, MIN]]
	correct_M = [[0, MIN, MIN, MIN, MIN, MIN], [MIN, MIN, MIN, MIN, MIN, MIN], [MIN, MIN, MIN, MIN, MIN, MIN], [MIN, MIN, MIN, MIN, MIN, MIN]]
	correct_v = [['end', 'left', 'left', 'left', 'left', 'left'], ['up', 'tbd', 'tbd', 'tbd', 'tbd', 'tbd'], ['up', 'tbd', 'tbd', 'tbd', 'tbd', 'tbd'], ['up', 'tbd', 'tbd', 'tbd', 'tbd', 'tbd']]

	class AmandaTmp:
		def __init__(self):
			self.gap_open = -3
			self.gap_ext = -1
			self.seqA = "ACACT" 
			self.seqB = "AAT"
			self.flag = "NW"
			self.sub_mat = "itdoesntmatter"

	[X, Y, M, val] = algs.distance_matrix(AmandaTmp(), method="testing_init")
	assert correct_X == X, "Failing to initialize X matrix properly"
	assert correct_Y == Y, "Failing to initialize Y matrix properly"
	assert correct_M == M, "Failing to initialize M matrix properly"
	assert correct_v == val, "Failing to initialize traceback matrix properly"

	#correct_X = [[-3, -inf, -inf, -inf, -inf, -inf], [-4, -inf, -inf, -inf, -inf, -inf], [-5, -2, -5, -6, -7, -8], [-6, -3, -3, -4, -7, -8]]
	#correct_Y = [[-3, -4, -5, -6, -7, -8], [-inf, -inf, -2, -3, -4, -5], [-inf, -inf, -5, -3, -4, -5], [-inf, -inf, -6, -6, -4, -5]]

	#[[[0, -inf, -inf, -inf, -inf, -inf], [-inf, 1, -2, -3, -4, -5], [-inf, -2, 0, -1, -4, -5], [-inf, -3, -3, -1, -2, -3]]]
	#M = algs.distance_matrix(AmandaTmp(), method="testing")
	#print(M)
	#assert False



def test_identical():
	nw = algs.NeedlemanWunsch(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM62")
	b = nw.align(algs.FastaRecord("sequences/prot-0004.fa"), algs.FastaRecord("sequences/prot-0004.fa"), print_flag=False)
	assert b.num_gaps == 0, "Failing base case"

	b = nw.align(algs.FastaRecord("sequences/prot-0008.fa"), algs.FastaRecord("sequences/prot-0008.fa"), print_flag=False)
	assert b.num_gaps == 0, "Failing base case"

	sw = algs.SmithWaterman(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM62")
	b = nw.align(algs.FastaRecord("sequences/prot-0004.fa"), algs.FastaRecord("sequences/prot-0004.fa"), print_flag=False)
	assert b.num_gaps == 0, "Failing base case"

	b = nw.align(algs.FastaRecord("sequences/prot-0008.fa"), algs.FastaRecord("sequences/prot-0008.fa"), print_flag=False)
	assert b.num_gaps == 0, "Failing base case"
	

def test_alignment_score():
	nw = algs.NeedlemanWunsch(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM50")
	b = nw.align(algs.FastaRecord("sequences/test1.fa"), algs.FastaRecord("sequences/test2.fa"), print_flag=False)
	assert b.num_gaps == 1, "Not finding gaps correctly"

	sw = algs.SmithWaterman(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM50")
	b = sw.align(algs.FastaRecord("sequences/test1.fa"), algs.FastaRecord("sequences/test2.fa"), print_flag=False)
	assert b.num_gaps == 1, "Not finding gaps correctly"

	sw = algs.SmithWaterman(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM62")
	b = sw.align(algs.FastaRecord("sequences/test5.fa"), algs.FastaRecord("sequences/test6.fa"), print_flag=False)
	assert b.top_score == 597, "Your scoring is messed up"
	assert b.top_pos == (112, 131), "Your position is messed up"






