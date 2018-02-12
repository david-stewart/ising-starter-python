from IsingLattice import IsingLattice
from ising_c import run_ising #import run_ising function from ising.py
import matplotlib.pyplot as plt

n=100
flip_prop=0.1

lattice = IsingLattice(n, flip_prop)

temp=2
num_steps=1250000
num_burnin= 50000
j=1
b=0.5

Msamp, Esamp = run_ising(lattice,temp,num_steps,num_burnin,j,b,disable_tqdm=False)

plt.figure(1)
plt.plot(Msamp)
plt.show()