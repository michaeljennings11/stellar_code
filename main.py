import numpy as np
from src.config import M
from src.initial_guess import get_guess
from src.integrate import delFit_xi, get_profiles
from scipy.optimize import minimize

print(f"M = {M} Msun")

print(f"starting solver-")
bnds = ((0, None), (0, None), (0, None), (0, 3))
found_vars = minimize(delFit_xi,x0=get_guess(M),bounds=bnds).x

profiles = get_profiles(found_vars)

np.savetxt("profiles.csv",profiles,delimiter=",",header="xi,F,p,r,T")
