#some functions for parametrized QNMs
import numpy as np
#import cmath as cm


def load_QNMs(field, what_l, what_n):
	"""
	loads Berti's QNMs for a given field, harmonic l and overtone n.
	prints an error message if files cannot be found.
	"""
	try:
		if field == "scalar":
			filename = "./QNM_GR/scalar/l" + str(what_l) + ".txt"
		elif field == "tensor_axial" or field == "tensor_polar":
			filename = "./QNM_GR/tensor/l" + str(what_l) + ".txt"
		w0_tmp = 2.0*np.loadtxt(filename, skiprows=0) #factor for units rH=1
		return w0_tmp[what_n]
	except:
		print("QNM file(s) not found, check if field="  + field + ", l=" + str(what_l) + ", n=" + str(what_n) + " is valid.")


def load_coefficients(field, what_l, what_n, what_precision):
	"""
	loads linear and quadratic coefficients for a given field, harmonic l and overtone n.
	prints an error message if files cannot be found.
	"""
	N_beta = 11 
	try:
		if field == "scalar":
			filename1 = "./QNM_coefficients/beyond_Schwarzschild/scalar/n" + str(what_n) + "_l" + str(what_l) + "_scalar_linear.txt"
			filename2 = "./QNM_coefficients/beyond_Schwarzschild/scalar/n" + str(what_n) + "_l" + str(what_l) + "_scalar_quadratic.txt"
		elif field == "tensor_axial":
			filename1 = "./QNM_coefficients/beyond_Schwarzschild/tensor_axial/n" + str(what_n) + "_l" + str(what_l) + "_axial_linear.txt"
			filename2 = "./QNM_coefficients/beyond_Schwarzschild/tensor_axial/n" + str(what_n) + "_l" + str(what_l) + "_axial_quadratic.txt"
		elif field == "tensor_polar":
			filename1 = "./QNM_coefficients/beyond_Schwarzschild/tensor_polar/n" + str(what_n) + "_l" + str(what_l) + "_polar_linear.txt"
			filename2 = "./QNM_coefficients/beyond_Schwarzschild/tensor_polar/n" + str(what_n) + "_l" + str(what_l) + "_polar_quadratic.txt"
		#linear
		d_       = np.zeros((N_beta, 2))
		T_tmp    = np.loadtxt(filename1, skiprows=0)
		for k in range(N_beta):
			idx = int(T_tmp[k][0])
			d_[idx] = np.array([T_tmp[k][1], T_tmp[k][2]])
		#quadratic
		e_       = np.zeros((N_beta, N_beta, 2))
		if what_precision == "quadratic":
			T_q_tmp  = np.loadtxt(filename2, skiprows=0)
			for k in range(len(T_q_tmp)):
				idx = int(T_q_tmp[k][0])
				idy = int(T_q_tmp[k][1])
				e_[idx, idy] = np.array([T_q_tmp[k][2], T_q_tmp[k][3]])
		return d_, e_
	except:
		print("Coefficient file(s) not found, check if field="  + field + ", l=" + str(what_l) + ", n=" + str(what_n) + " is valid.")


def load_coefficient_file(what_l, what_m, what_n, what_k):
    """
    loads linear coefficients for a given harmonic l, harmonic m, overtone n, and order k
    prints an error message if file cannot be found.
    """
    try:
        if what_k < 0 and what_m < 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "_l" + str(what_l) + "_mm" + str(-what_m) + "_km" + str(-what_k) + ".txt"
        if what_k < 0 and what_m >= 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "_l" + str(what_l) + "_m" + str(what_m) + "_km" + str(-what_k) + ".txt"
        if what_k >= 0 and what_m < 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "_l" + str(what_l) + "_mm" + str(-what_m) + "_k" + str(what_k) + ".txt"
        if what_k >= 0 and what_m >= 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "_l" + str(what_l) + "_m" + str(what_m) + "_k" + str(what_k) + ".txt"
        #print(filename)
        T_tmp     = np.genfromtxt(filename)
        return T_tmp
    except:
    	print("Coefficient file(s) not found, check if l=" + str(what_l) + ", m=" + str(what_m) + ", n=" + str(what_n) + " is valid." + "\n" + filename)


def load_diag_coefficient_file(what_l, what_m, what_n, what_k):
    """
    loads quadratic diagonal coefficients for a given harmonic l, harmonic m, overtone n, and order k
    prints an error message if file cannot be found.
    """
    try:
        if what_k < 0 and what_m < 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "l" + str(what_l) + "mm" + str(-what_m) + "km" + str(-what_k) + "_quadratic.txt"
        if what_k < 0 and what_m >= 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "l" + str(what_l) + "m" + str(what_m) + "km" + str(-what_k) + "_quadratic.txt"
        if what_k >= 0 and what_m < 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "l" + str(what_l) + "mm" + str(-what_m) + "k" + str(what_k) + "_quadratic.txt"
        if what_k >= 0 and what_m >= 0:
            filename = "./QNM_coefficients/beyond_Teukolsky/l" + str(what_l) + "/n" + str(what_n) + "l" + str(what_l) + "m" + str(what_m) + "k" + str(what_k) + "_quadratic.txt"
        #print(filename)
        T_tmp     = np.genfromtxt(filename)
        return T_tmp
    except:
    	print("Coefficient file(s) not found, check if l=" + str(what_l) + ", m=" + str(what_m) + ", n=" + str(what_n) + " is valid." + "\n" + filename)
    	filename = "./QNM_coefficients/beyond_Teukolsky/l2/n0l2m2k0_quadratic.txt"
    	T_tmp    = np.genfromtxt(filename)
    	return (1e-15)*T_tmp
