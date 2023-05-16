from src.config import M
from src.initial_guess import get_guess
from src.integrate import delFit_xi
from scipy.optimize import minimize

print(f"M = {M} Msun")

print(f"starting solver-")
bnds = ((0, None), (0, None), (0, None), (0, 3))
minimize(delFit_xi,x0=get_guess(M),bounds=bnds)

