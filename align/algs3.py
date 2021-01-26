from align import algs
MIN=-float("inf")
correct_X = [[-3, MIN, MIN, MIN, MIN, MIN], [-4, MIN, MIN, MIN, MIN, MIN], [-5, -3, -9, -8, -11, -12], [-6, -4, -4, -6, -9, -10]]
correct_Y = [[-3, -4, -5, -6, -7, -8], [MIN, MIN, -3, -4, -5, -6], [MIN, MIN, -7, -4, -5, -6], [MIN, MIN, -10, -8, -5, -6]]
correct_M = [[0, MIN, MIN, MIN, MIN, MIN], [MIN, 1, -5, -4, -7, -8], [MIN, -3, 0, -2, -5, -6], [MIN, -6, -4, -1, -3, -4]]
	
class AmandaTmp:
	def __init__(self):
		self.gap_open = -3
		self.gap_ext = -1
		self.seqA = "AAT"
		self.seqB = "ACACT"
		self.flag = "NW"
		self.sub_mat = "itdoesntmatter"

a = AmandaTmp()
print(a.gap_open)

#[X, Y] = algs.distance_matrix(AmandaTmp(), method="testing")