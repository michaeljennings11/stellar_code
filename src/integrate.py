import numpy as np
from scipy.integrate import odeint
from .config import star_params,xi_delta_core
from .boundaries import load1,load2
from .diff_eqs import derivs_xi

xi_vals1 = np.logspace(np.log10(xi_delta_core),np.log10(0.999),5_000)
xi_vals2 = 1.0 - np.logspace(-9,-3,1000)

def delFit_xi(guesses):
    pc_guess,Tc_guess,Fs_guess,Rs_guess = guesses
    
    Fc,pc,rc,Tc = load1(pc_guess,Tc_guess,star_params,xi_delta_core)
    Fs,ps,rs,Ts = load2(Rs_guess,Fs_guess,star_params)
    
    res1 = odeint(func=derivs_xi, y0=[Fc,pc,rc,Tc], t=xi_vals1)
    res2 = odeint(func=derivs_xi, y0=[Fs,ps,rs,Ts], t= xi_vals2)
    score = np.abs(np.subtract(res1[-1],res2[-1]))
    return np.sum(score)
