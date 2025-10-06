# stellar_code

stellar_code is a 1D stellar structure modelling code written for the final project
of the Stellar Evolution graduate level course at Johns Hopkins University.

Given a star's mass, hydrogen, helium and metal mass fraction,
stellar_code will solve the set of 4 differential equations for the radial profiles.

The full instructions for the final project given in the course are found in instructions.pdf.

## run

Use the following to run:

```shell
# run stellar code
python main.py

# main.py will produce csv file called profiles.csv with xi,F,p,r,T variables

# change stellar mass and composition in config.py
```
