import numpy as np
from .constants import N_A, k, a

def mu_ionized(X):
	return 4/(3+5*X)

def q_eos(p,T,mu):
	return np.sqrt((mu/(N_A*k*T))*(p**4 - (a/3.0)*T**4))

def Pgas(T,rho,mu):
	return rho*N_A*k*T/mu

def Prad(T):
	return a*T**4/3

def P_total(T,rho,mu):
	return Pgas(T,rho,mu) + Prad(T)
