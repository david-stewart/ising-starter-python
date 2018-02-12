from IsingLattice import IsingLattice
from ising_c import run_ising #import run_ising function from ising.py
from ising_old import run_ising_old
import matplotlib.pyplot as plt

n=50
flip_prop=0.1

lattice = IsingLattice(n, flip_prop)

temp=2
num_steps=100000
num_burnin= 50000
j=1
b=0

Msamp, Esamp = run_ising(lattice,temp,num_steps,num_burnin,j,b,disable_tqdm=False)

num_steps=100000

Msamp_old, Esamp_old, spin = run_ising_old(n,temp,num_steps,num_burnin,flip_prop,j,b,disable_tqdm=False)

plt.figure(1)
plt.plot(Msamp)
plt.figure(2)
plt.plot(Msamp_old)
plt.show()