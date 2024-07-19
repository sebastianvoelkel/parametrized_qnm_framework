#file that contains the QNM_class
import numpy as np
import cmath as cm
import extra_functions as EF
from scipy.interpolate import CubicSpline

class QNM_class:
    """
    This class contains all the functions to compute QNMs and separation constants.
    GR QNMs for Kerr are read and interpolated from file taken from Berti's ringdown website.
    Linear coefficients are read from file and interpolated with cubic splines.

    Computation of QNMs and separation constants:
    #call class
    QNM = QNM_class()
    #initialize field content and approximant
    QNM.initialize("tensor", "linear")
    #compute QNMs for X (array length = 11)
    QNM.qnm(l,m,n,a,X)
    """

    
    def __init__(self):
        self.omega_GR     = None
        self.B_GR         = None
        self.d_omega      = None
        self.d_B          = None
        self.e_omega_diag = None
        self.e_B_diag     = None
        self.field        = None
        self.approximant  = None

    
    def set_field(self, field):
        if field in ["tensor", "scalar", "tensor_axial", "tensor_polar"]: 
            self.field = field
        else:
            print(f"field = {field} is not a valid option.")

            
    def set_approximant(self, approximant):
        if approximant in ["linear", "quadratic"]: 
            self.approximant = approximant
        else:
            print(f"approximant = {approximant} is not a valid option.")
   

    def build_splines_GR(self):
        """
        load GR QNMs and separation constant for Kerr
        compute cubic splines for given [l,m,n,case]
        """
        omega_GR_spline_dict = {}
        B_GR_spline_dict     = {}
        for l in [2,3,4]:
            for m in range(-l,l+1):
                for n in [0,1,2]:
                    #what_n +1 because files n=0 is n=1
                    if m < 0:
                        filename = "./QNM_GR/tensor_Kerr/l" + str(l) + "/n" + str(n+1) + "l" + str(l) + "mm" + str(-m) + ".dat"
                    if m >= 0:
                        filename = "./QNM_GR/tensor_Kerr/l" + str(l) + "/n" + str(n+1) + "l" + str(l) + "m" + str(m) + ".dat"
                    T_tmp = np.loadtxt(filename, skiprows=0)
                    omega_GR_spline_dict[l, m, n, 0] = CubicSpline(T_tmp[:,0]/2., 2.*T_tmp[:,1])
                    omega_GR_spline_dict[l, m, n, 1] = CubicSpline(T_tmp[:,0]/2., 2.*T_tmp[:,2])
                    B_GR_spline_dict[l, m, n, 0]     = CubicSpline(T_tmp[:,0]/2., T_tmp[:,3])
                    B_GR_spline_dict[l, m, n, 1]     = CubicSpline(T_tmp[:,0]/2., T_tmp[:,4])
        self.omega_GR = omega_GR_spline_dict
        self.B_GR     = B_GR_spline_dict

    
    def build_splines_coefficients(self):
        """
        load linear coefficient table
        compute cubic splines for given [l,m,n,case,k]
        """
        d_omega_spline_dict        = {}
        d_B_spline_dict            = {}
        e_omega_spline_diag_dict   = {}
        e_B_spline_diag_dict       = {}
        for l in [2,3,4]:
            for m in range(-l,l+1):
                for n in [0,1,2]:
                    for case in [0,1]:#real, imaginary
                        for k in range(-6,5):
                            coeffs      = EF.load_coefficient_file(l, m, n, k)
                            coeffs_diag = EF.load_diag_coefficient_file(l, m, n, k)
                            x      = coeffs[~np.isnan(coeffs).any(axis=1)] #delete all rows with nans
                            x_diag = coeffs_diag[~np.isnan(coeffs_diag).any(axis=1)] 
                            d_omega_spline_dict[l, m, n, k, case]      = CubicSpline(x[:,0], x[:,case+1])#, bc_type="natural"
                            d_B_spline_dict[l, m, n, k, case]          = CubicSpline(x[:,0], x[:,case+3])#, bc_type="natural"
                            e_omega_spline_diag_dict[l, m, n, k, case] = CubicSpline(x_diag[:,0], x_diag[:,case+1])#, bc_type="natural"
                            e_B_spline_diag_dict[l, m, n, k, case]     = CubicSpline(x_diag[:,0], x_diag[:,case+3])#, bc_type="natural"
        self.d_omega      = d_omega_spline_dict
        self.d_B          = d_B_spline_dict
        self.e_omega_diag = e_omega_spline_diag_dict
        self.e_B_diag     = e_B_spline_diag_dict

    
    def qnm_rot(self, what_l, what_m, what_n, what_a, X):
    	"""
    	computes linear QNMs for rotating BH for given set of X = [alpha_k]
    	"""
    	w0 = np.array([self.omega_GR[what_l, what_m, what_n, 0](what_a), self.omega_GR[what_l, what_m, what_n, 1](what_a)])
    	wr_res, wi_res = 0., 0.
    	index_shift    = -6
    	for i in range(len(X)):
            d_temp_i = self.d_omega[what_l, what_m, what_n, i+index_shift, 0](what_a) + 1j*self.d_omega[what_l, what_m, what_n, i+index_shift, 1](what_a)
            X_temp_i = X[i].real + 1j*X[i].imag
            wr_res   = wr_res + (d_temp_i*X_temp_i).real
            wi_res   = wi_res + (d_temp_i*X_temp_i).imag
    	w_res = w0 + np.array([wr_res, wi_res])
    	return w_res


    def qnm_rot_quad_diag(self, what_l, what_m, what_n, what_a, X):
        """
        computes "diagonal" quadratic QNM corrections for given set of X = [alpha_k]
        """
        wr_res, wi_res = 0., 0.
        index_shift    = -6
        for i in range(len(X)):
    	    for j in range(len(X)):
                #attention for summing twice i,j, make i=j case, load diag coeff
                if i==j:
                    e_temp_ij = self.e_omega_diag[what_l, what_m, what_n, i+index_shift, 0](what_a) + 1j*self.e_omega_diag[what_l, what_m, what_n, i+index_shift, 1](what_a)
                else:
                    e_temp_ij = 0
                X_temp_i  = X[i].real   + 1j*X[i].imag
                X_temp_j  = X[j].real   + 1j*X[j].imag
                wr_res    = wr_res      + 0.5*(e_temp_ij*X_temp_i*X_temp_j).real
                wi_res    = wi_res      + 0.5*(e_temp_ij*X_temp_i*X_temp_j).imag
        w_res = np.array([wr_res, wi_res])
        w_lin = self.qnm_rot(what_l, what_m, what_n, what_a, X)
        return w_res + w_lin


    def B_rot(self, what_l, what_m, what_n, what_a, X):
    	"""
    	computes linear separation constant for rotating BH for given set of X = [alpha_k]
    	"""
    	B0 = np.array([self.B_GR[what_l, what_m, what_n, 0](what_a), self.B_GR[what_l, what_m, what_n, 1](what_a)])
    	Br_res, Bi_res = 0., 0.
    	index_shift    = -6
    	for i in range(len(X)):
            d_temp_i = self.d_B[what_l, what_m, what_n, i+index_shift, 0](what_a) + 1j*self.d_B[what_l, what_m, what_n, i+index_shift, 1](what_a)
            X_temp_i = X[i].real + 1j*X[i].imag
            Br_res   = Br_res + (d_temp_i*X_temp_i).real
            Bi_res   = Bi_res + (d_temp_i*X_temp_i).imag
    	B_res = B0 + np.array([Br_res, Bi_res])
    	return B_res


    def compute_qnm_lin(self, field, what_l, what_n, X):
    	"""
    	computes linear QNM corrections for given set of X = [alpha_k]
    	"""
    	w0     = EF.load_QNMs(field, what_l, what_n)
    	d_, e_ = EF.load_coefficients(field, what_l, what_n, what_precision="linear")
    	wr_res, wi_res = 0., 0.
    	for i in range(len(X)):
    		d_temp_i = d_[i][0]  + 1j*d_[i][1]
    		X_temp_i = X[i].real + 1j*X[i].imag
    		wr_res   = wr_res + (d_temp_i*X_temp_i).real
    		wi_res   = wi_res + (d_temp_i*X_temp_i).imag
    	w_res = w0 + np.array([wr_res, wi_res])
    	return w_res


    def compute_qnm_quad(self, field, what_l, what_n, X):
    	"""
    	computes quadratic QNM corrections for given set of X = [alpha_k]
    	"""
    	w0     = EF.load_QNMs(field, what_l, what_n)
    	d_, e_ = EF.load_coefficients(field, what_l, what_n, what_precision="quadratic")
    	wr_res, wi_res = 0., 0.
    	for i in range(len(X)):
    		for j in range(len(X)):
    			e_temp_ij = e_[i][j][0] + 1j*e_[i][j][1]
    			X_temp_i  = X[i].real   + 1j*X[i].imag
    			X_temp_j  = X[j].real   + 1j*X[j].imag
    			wr_res    = wr_res      + 0.5*(e_temp_ij*X_temp_i*X_temp_j).real
    			wi_res    = wi_res      + 0.5*(e_temp_ij*X_temp_i*X_temp_j).imag
    	w_res = np.array([wr_res, wi_res])
    	w_lin = self.compute_qnm_lin(field, what_l, what_n, X)
    	return w_res + w_lin
    
    
    def qnm_nonrot(self, field, what_l, what_n, X, what_precision):
    	"""
    	computes QNM corrections (linear or quadratic) for given set of X = [alpha_k]
    	"""
    	if what_precision == 'linear':
    		res = self.compute_qnm_lin(field, what_l, what_n, X)
    	if what_precision == 'quadratic':
    		res = self.compute_qnm_quad(field, what_l, what_n, X)
    	return res

    
    def qnm(self, what_l, what_m, what_n, what_a, X):
        """
        wrapper for computing QNMs (Schwarzschild or Kerr)
        """
        if self.field == "tensor":
            return self.qnm_rot(what_l, what_m, what_n, what_a, X)
        if self.field in ["scalar", "tensor_axial", "tensor_polar"]:
            return self.qnm_nonrot(self.field, what_l, what_n, X, self.approximant)  


    def qnm_rot_error(self, what_l, what_m, what_n, what_a, X):
        """
        wrapper for returning error of quadratic-linear qnm
        """
        if self.field == "tensor":
            w_linear    = self.qnm_rot(what_l, what_m, what_n, what_a, X)
            w_quad_diag = self.qnm_rot_quad_diag(what_l, what_m, what_n, what_a, X)
            #print(w_linear, w_quad_diag)
            return w_quad_diag - w_linear

    
    def constant(self, what_l, what_m, what_n, what_a, X):
        """
        wrapper for computing separation constant (only Kerr)
        """
        if self.field == "tensor":
            return self.B_rot(what_l, what_m, what_n, what_a, X)


    def initialize(self, field, approximant):
        self.set_field(field)
        self.set_approximant(approximant)
        if approximant == "linear" and field == "tensor":
            self.build_splines_GR()
            #self.build_splines_omega_GR()
            self.build_splines_coefficients()
            print(f"Initialized QNM class for computations with field = {field} and approximant = {approximant}.")
        elif approximant in ["linear", "quadratic"] and field in ["scalar", "tensor_axial", "tensor_polar"]:
            print("Nonrotating case.")
        else:
            print(f"Initialization failed, check field = {field} and approximant = {approximant}")
