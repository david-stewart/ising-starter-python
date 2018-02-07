import random
import numpy as np
from scipy import signal
from annealing import B_anneal, T_anneal
from ctypes import *

try:
    __IPYTHON__
except:
    from tqdm import tqdm

def run_ising(c_matrix, N,T,num_steps,num_burnin,flip_prop,J,B,disable_tqdm=False):
    # c_matrix.print_spins()
    c_matrix.set_flip_prop(c_float(flip_prop))
    c_matrix.set_J(c_float(J))
    c_matrix.rand_spins()
    # Description of parameters:
    # N = Grid Size
    # T = Temperature (normalized to k_B = 1)
    # num_steps = Number of steps to run in total (including burnin steps)
    # num_burnin = Number of steps to use for the burnin process.  This isn't
    #	used in this code but you might need it if you change this to try and
    #	get better convergence.
    # J = Interaction strength
    # B = Applied magnetic field
    # flip_prop = Total ratio of spins to possibly flip per step

    Msamp, Esamp = [],[] #Arrays to hold magnetization and energy values

    # We obtain the sum of nearest neighbors by convoluting
    #   this matrix with the spin matrix
    try:
        __IPYTHON__
        steps = range(num_steps)
    except:
        if disable_tqdm:
            steps = range(num_steps)
        else:
            steps = tqdm(range(num_steps))

    # Evolve the system
    for step in steps:

        try:
            __IPYTHON__
        except:
            if disable_tqdm:
                pass
            else:
                steps.set_description("Working on T = %.2f" % T)
                steps.refresh() # to show immediately the update

        #implement annealing in annealing.py file
        T_step = T_anneal(T, step, num_steps, num_burnin)
        B_step = B_anneal(B, step, num_steps, num_burnin)

        c_matrix.step(c_float(T_step), c_float(B_step))
        # c_matrix.print_E()
        Msamp.append(float(c_matrix.get_M()))
        Esamp.append(float(c_matrix.get_E()))

    return Msamp, Esamp
