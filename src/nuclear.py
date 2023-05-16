import numpy as np

def XCNO(Z):
	return 0.70*Z

def zeta(X,Y):
    return 2*X+(3/2)*Y

def f_weak(T,rho,X,Y,Z1=1,Z2=1):
    T7 = T/1e7
    return np.exp(5.92e-3*Z1*Z2*(zeta(X,Y)*rho/(T7**3))**(1/2))

def eps_pp(T,rho,X,Y,psi=1):
    T9 = T/1e9
    g11 = 1.+3.82*T9+1.51*T9**2+0.144*T9**3-0.0114*T9**4
    return 2.57e4*psi*f_weak(T,rho,X,Y)*g11*rho*(X**2)*(T9**(-2/3))*np.exp(-3.381/(T9**(1/3)))

def eps_CNO(T,rho,X,Z):
    T9 = T/1e9
    g14_1 = (1.-2*T9+3.41*T9**2-2.43*T9**3)
    return 8.24e25*g14_1*rho*XCNO(Z)*X*(T9**(-2/3))*np.exp(-15.231*(T9)**(-1/3)-(T9/0.8)**2)

def eps_total(T,rho,X,Y,Z,psi=1):
    ept = eps_pp(T,rho,X,Y,psi=1) + eps_CNO(T,rho,X,Z)
    return ept


