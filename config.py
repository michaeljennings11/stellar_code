from .opacity import interp_opac
from .eos import mu_ionized

opal_file = '../opacities/GN93/GN93hz'

M = 3
X = 0.7
Y = 0.28
Z = 1-X-Y
mu = mu_ionized(X)
star_params = (M,X,Y,Z)

interp = interp_opac((X,Y,0,0),opal_file,method="regular")
