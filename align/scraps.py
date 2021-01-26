def distance_matrix(seqA, seqB, gap_open, gap_ext, sub_mat, method="real"):
    MIN=-float("inf")
    dim_i = len(seqA) + 1
    dim_j = len(seqB) + 1

    X = [[_init_y(i, j, gap_open, gap_ext) for i in range(0, dim_i)] for j in range(0, dim_j)]
    Y = [[_init_x(i, j, gap_open, gap_ext) for i in range(0, dim_i)] for j in range(0, dim_j)]
    M = [[_init_m(i, j) for i in range(0, dim_i)] for j in range(0, dim_j)]
    val = [[_init_v(i, j) for i in range(0, dim_i)] for j in range(0, dim_j)]

    max_val = 0
    max_pos = None
    #print(M); print(X); print(Y); print(val)

    for i in range(1, dim_j):
        for j in range(1, dim_i):
            tmp = _match(sub_mat, seqA[j - 1], seqB[i - 1], method) + M[i-1][j-1]
            Y[i][j] = max(
                          (M[i][j-1] + gap_open),
                          (Y[i][j-1] + gap_ext))
            X[i][j] = max(
                          (M[i-1][j] + gap_open), 
                          (X[i-1][j] + gap_ext)
                          )

            #ugh = [tmp, X[i][j], Y[i][j]]
            ugh = [tmp, X[i][j], Y[i][j], 0]
            M[i][j]  = max(ugh)

            if M[i][j] > max_val:                 #calculate it so that we know where to start traceback
                max_val = M[i][j]
                max_pos = (i,j)

            maxmat = ugh.index(max(ugh))
            if maxmat == 0:
                val[i][j] = "diag"
            elif maxmat == 1:
                val[i][j] = "up"
            elif maxmat == 2: 
                val[i][j] = "left"
            elif maxmat == 3:
                val[i][j] = "end"
            else:
                print("ugh bro howd this happen"); exit()
    #print(M); print(X); print(Y); print(val)
    return [val, M, max_val, max_pos]
    


def backtrace(seqA, seqB, X, Y, M, sub_mat, method="real"):
    sequ1 = ''
    sequ2 = ''

    #intialize at last corner
    i = len(seqB); j = len(seqA)

    while (i>0 or j>0):

        ugh = [ M[i][j], X[i][j], Y[i][j]]
        #print("i=", i); print("j=", j); print("M= ", M[i][j]) #; print("Mi-1,j-1 = ", M[i-1][j-1])
        #print("X =", X[i][j]); print("Y =", Y[i][j])
        #print("ugh= ", ugh)
        maxmat = ugh.index(max(ugh))
        #print(maxmat)
        if maxmat == 0 and i>0 and j>0: #diag
            #note, ugh.index(max(ugh)) will return the first occurence first if tie, so important to keep M[i,j] first in list
            sequ1 += seqA[j-1]
            sequ2 += seqB[i-1]
            i -= 1; j -= 1
        elif maxmat == 1 and i > 0: #down
            sequ1 += '-'
            sequ2 += seqB[i-1]
            i -= 1
        elif maxmat == 2 and j>0: #up
            sequ1 += seqA[j-1]
            sequ2 += '-'
            j -= 1
        else:
            print("How tf did we get here")
            print(maxmat); print(i, j)
            break
            #exit(1)

    if len(sequ1) != len(sequ2):
        raise ValueError('issue in traceback: alignments not the same length')

    a = ''.join(reversed(sequ1)), ''.join(reversed(sequ2))
    return a




def backtrace(seqA, seqB, X, Y, M, sub_mat, method="real"):
    sequ1 = ''
    sequ2 = ''

    #intialize at last corner
    i = len(seqB); j = len(seqA)

    while (i>0 or j>0):
        ugh = [ M[i][j], X[i][j], Y[i][j]]
        maxmat = ugh.index(max(ugh))
        if maxmat == 0:
            #note, ugh.index(max(ugh)) will return the first occurence first if tie, so important to keep M[i,j] first in list
            #diag
            #print("diag")
            sequ1 += seqA[j-1]
            sequ2 += seqB[i-1]
            i -= 1; j -= 1
        elif maxmat == 1:
            #print("up")
            sequ1 += seqA[j-1]
            sequ2 += '-'
            j -= 1
        elif maxmat == 2:
            #print("down")
            sequ1 += '-'
            sequ2 += seqB[i-1]
            i -= 1
        else:
            print("How tf did we get here")
            exit(1)

    a = ''.join(reversed(sequ1)), ''.join(reversed(sequ2))
    return a    



def backtrace(seqA, seqB, X, Y, M, sub_mat, method="real"):
    sequ1 = ''
    sequ2 = ''

    #intialize at last corner
    i = len(seqB); j = len(seqA)
    
    while ((i>0 or j>0)):
        ugh = max([ M[i][j], X[i][j], Y[i][j]])
        if (i>0 and j>0 and ugh == M[i-1][j-1] + _match(sub_mat, seqA[j - 1], seqB[i - 1], method)):
            sequ1 += seqA[j-1]
            sequ2 += seqB[i-1]
            i -= 1; j -= 1
        elif (i>0 and ugh == X[i][j]):
            sequ1 += '_'
            sequ2 += seqB[i-1]
            i -= 1
        elif (j>0 and ugh == Y[i][j]):
            sequ1 += seqA[j-1]
            sequ2 += '_'
            j -= 1
        else:
            print("geezus fuck")
            print(i)
            print(j)
            break

    #a = ''.join(reversed(sequ1)), ''.join(reversed(sequ2))
    sequ1r = ' '.join([sequ1[j] for j in range(-1, -(len(sequ1)+1), -1)])
    sequ2r = ' '.join([sequ2[j] for j in range(-1, -(len(sequ2)+1), -1)])
    print(sequ1r)
    print(sequ2r)
    return sequ1r, sequ2r





        #diag = M[i-1][j-1]
        #print("i=", i)
        #print("j=", j)
        #print(M[i][j])
        #print(M[i-1][j-1])
        #print(_match(sub_mat, seqA[j-1], seqB[i-1], method))
        #MOVE DIAGNOL
        if (i>0 and j>0 and cell == diag + _match(sub_mat, seqA[j-1], seqB[i-1], method)):
            print("Diag")
            sequ1 += seqA[j-1]
            sequ2 += seqB[i-1]
            i -= 1; j -= 1
        #MOVE UP # subtract from columns, but not row. 
        elif (i>0 and cell == Y[i][j]):
            print("up")
            sequ1 += seqA[j-1]
            sequ2 += '-'
            i -= 1
        #MOVE DOWN # subtract from columns, but not row.
        elif (j>0 and cell == X[i][j]):
            print("down")
            sequ1 += '-'
            sequ2 += seqB[i-1]
            j -= 1
        else:
            print("How tf did we get here")
            print(i,j,cell, diag, Y[i][j],X[i][j])
            exit(1)

#class SmithWaterman(PairwiseAligner):
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
	#pass

#class NeedlemanWunsch(PairwiseAligner):
	#https://anhaidgroup.github.io/py_stringmatching/v0.3.x/_modules/py_stringmatching/similarity_measure/needleman_wunsch.html#NeedlemanWunsch
	#pass



def alignment_string(aligned_seq1, aligned_seq2):
    '''Construct a special string showing identities, gaps, and mismatches.
    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.
    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    '''
    # Build the string as a list of characters to avoid costly string
    # concatenation.
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

    return ''.join(alignment_string), idents, gaps, mismatches





    def read_fasta(sequencefile):
    #OneEntryCheck = T
    name, seq = None, []
    #with open('../' + sequencefile) as fp:
    for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                print(name)
                #OneEntryCheck = T
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))
            

        #final_seq = ''.join(seq)
        #if final_seq.isalpha():
        #   final_seq = final_seq.upper()
        #else:
        #   print("You've entered an invalid sequence. Please try again.")
        #   exit()
    #return name, final_seq




