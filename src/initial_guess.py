import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from .config import M,mu
from .boundaries import p_central,T_central
from .constants import Lsun

def RTL_lookup(M):
    stars_df = pd.read_csv('./EEM_dwarf_UBVIJHK_colors_Teff.txt',skiprows=22,nrows=118,delim_whitespace=True)
    stars_df = stars_df.replace('...',float("nan"))
    interp_MR = interp1d(stars_df["Msun"].astype(float),stars_df["R_Rsun"].astype(float))
    interp_MTeff = interp1d(stars_df["Msun"].astype(float),stars_df["Teff"].astype(float))
    interp_ML = interp1d(stars_df["Msun"].astype(float),stars_df["logL"].astype(float))
    return interp_MR(M),interp_MTeff(M),np.power(10,interp_ML(M))

def get_guess(M):
		Rs_guess,_,L_guess = RTL_lookup(M)
		Fs_guess = L_guess*Lsun
		pc_guess = p_central(M,Rs_guess)
		Tc_guess = T_central(M,Rs_guess,mu)
		return [pc_guess,Tc_guess,Fs_guess,float(Rs_guess)]
