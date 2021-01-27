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
		def __init__(self, sequence):
			self.name = "adlkf"
			self.seq = sequence
	nw = algs.NeedlemanWunsch(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM62")
	a= nw.align(AmandaTmp("ACACT"), AmandaTmp("AAT"), print_flag=False)
	b = a.initialize_matrices()

	#[X, Y, M, val] = algs.distance_matrix(AmandaTmp(), method="testing_init")
	assert correct_X == b.X, "Failing to initialize X matrix properly"
	assert correct_Y == b.Y, "Failing to initialize Y matrix properly"
	assert correct_M == b.M, "Failing to initialize M matrix properly"
	assert correct_v == b.TB, "Failing to initialize traceback matrix properly"



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
	class AmandaTmp:
		def __init__(self, sequence):
			self.name = "adlkf"
			self.seq = sequence

	nw = algs.NeedlemanWunsch(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM62")
	sw = algs.SmithWaterman(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM62")
	b = nw.align(AmandaTmp("YYYWWW"), AmandaTmp("YYYAAA"), print_flag=False)
	assert b.num_matches == 3, "Not finding matches correctly -- NW"
	b = sw.align(AmandaTmp("YYYWWW"), AmandaTmp("YYYAAA"), print_flag=False)
	assert b.num_matches == 3, "Not finding matches correctly -- SW"

	b = nw.align(AmandaTmp("YYYWWWAAA"), AmandaTmp("YYYAAA"), print_flag=False)
	assert b.num_gaps == 3, "Not finding gaps correctly -- NW"
	b = sw.align(AmandaTmp("YYYWWWAAA"), AmandaTmp("YYYAAA"), print_flag=False)
	assert b.num_gaps == 3, "Not finding gaps correctly --SW"

	nw = algs.NeedlemanWunsch(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM50")
	sw = algs.SmithWaterman(gap_open=-3, gap_extension = -1, substitutionMatrix = "BLOSUM50")
	#b = nw.align(AmandaTmp("CALM"), AmandaTmp("ACALMA"), print_flag=True)
	#print(b.top_score)

	#assert b.top_score == 22, "Not generating scores correctly -- NW"
	b = sw.align(AmandaTmp("CALM"), AmandaTmp("ACALMA"), print_flag=False)
	assert b.top_score == 30, "Not generating scores correctly -- NW"

	b = nw.align(AmandaTmp("ACACT"), AmandaTmp("AAT"), print_flag=False)
	assert b.top_score == 7, "Not generating scores correctly -- NW"
	b = sw.align(AmandaTmp("ACACT"), AmandaTmp("AAT"), print_flag=False)
	assert b.top_score == 9, "Not generating scores correctly -- NW"

