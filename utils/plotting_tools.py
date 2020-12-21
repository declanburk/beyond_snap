import csv
import numpy as np

def poly_at_time(p, t):
	return np.sum(np.array([p[i] * (t ** i) for i in range(len(p))]))

def deriv_coeff(n, d):
    if (d == 0):
        return np.ones(n)
    else:
        ind = np.array([max(0, i - d + 1) for i in range(n)])
        rem = deriv_coeff(n, d - 1)
        return np.multiply(ind, rem)

def deriv_at_time(p, d, t):
	tim = np.array([t ** (i - d) if i - d >= 0 else 0 for i in range(n)])
	bas = deriv_coeff(n, d)
	return np.dot(p, np.multiply(tim, bas))

def eval_poly_coeffs(p, t, k, n, m):
	T = []
	P = []
	for i in range(k):
		TT = [m * j for j in range(int(1 / m))]
		for tt in TT:
			temp = tt * (t[i+1] - t[i]) + t[i]
			T.append(temp)
			P.append(poly_at_time(p[i * n : (i + 1) * n], temp))
	return T, P

def eval_normalised_poly_coeffs(p, t, k, n, m):
	T = []
	P = []
	for i in range(k):
		TT = [m * j for j in range(int(2 / m))]
		for tt in TT:
			T.append(((t[i + 1] - t[i]) / 2) * (tt - 1) + ((t[i + 1] + t[i]) / 2))
			P.append(poly_at_time(p[i * n : (i + 1) * n], tt - 1))
	return T, P

def read_spline_from_csv(name):
	t = []; p = [];	b = []
	with open(name, 'rt') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=' ', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
		for ir, row in enumerate(csvreader):
			if (ir == 0):
				k = int(row[0])
			elif (ir == 1):
				for r in row:
					if (r != ''):
						t.append(r)
			elif (ir == 2):
				for r in row:
					if (r != ''):
						b.append(r)
			elif (ir == 3):
				J = float(row[0])
			else:
				for r in row:
					if (r != ''):
						p.append(r)
		return p, t, b, k, J