# Python script to solve the PDE model for biomass and age classes
# Created by Luca Rossini
# e-mail: luca.rossini@unitus.it
# Last update: 25 July 2021


# List of import

import pandas as pd
import numpy as np
import sys
import random as rd
from random import seed
from random import random
from Parameters import *


# Generating random numbers:

# Seeds of the random generator

seed = rd.randrange(sys.maxsize) + N_iter
rd.seed(sys.maxsize)

# Briere

def BriPar_generator(a, T_L, T_M, m):
    a = a + (random() * 2 * var * a)
    T_L = T_L + (random() * 2 * var * T_L )
    T_M = T_M + (random() * 2 * var * T_M )
    m = m + (random() * 2 * var * m )

    BriPar = [a, T_L, T_M, m]
    return BriPar


# Biomass

def BioPar_generator(a_1, h_b, P, var): 
    a_1 = a_1 + (random() * 2 * var * a_1 )
    h_b = h_b + (random() * 2 * var * h_b )
    P = P + (random() * 2 * var * P )

    BioPar = [a_1, h_b, P]
    return BioPar


# Fertility

def FertPar_generator(gamma, sigma_T, T_star, T_low, T_max, var):
    gamma = gamma + (random() * 2 * var * gamma )
    sigma_T = sigma_T + (random() * 2 * var * sigma_T )
    T_star = T_star + (random() * 2 * var * T_star )
    T_low = T_low + (random() * 2 * var * T_low )
    T_max = T_max + (random() * 2 * var * T_max )

    FertPar = [gamma, sigma_T, T_star, T_low, T_max]
    return FertPar


# Mortality

def MortPar_generator(a_2, T_opt, var):
    a_2 = a_2 + (random() * 2 * var * a_2 )
    T_opt = T_opt + (random() * 2 * var * T_opt )

    MortPar = [a_2, T_opt]
    return MortPar


# Definition of the temperature function

def TempFunc(TemPar, i):
    return TemPar[0] + TemPar[1] * np.cos(((2 * np.pi)/TemPar[2]) * i)


# Definition of Briere rate function

def BriFunc(BriPar, TemPar, i):
    return BriPar[0] * TempFunc(TemPar, i) * (TempFunc(TemPar, i) - BriPar[1]) * pow((BriPar[2] - TempFunc(TemPar, i)), (1 / BriPar[3]))


# Definition of the functional response

def BioFunc(BioPar):
    return (BioPar[0] * BioPar[2]) / (1 + BioPar[0] * BioPar[1] * BioPar[2])


# Definition of the fertility rate function

def FertFunc(FertPar, TemPar, i):
    FP_1 = FertPar[0] * BriFunc(BriPar, TemPar, i)
    FP_2 = 1 / (pow( 2 * np.pi, 0.5) * FertPar[1])
    FP_3 = np.exp( -(  pow(TempFunc(TemPar, i) - FertPar[2], 2) ) / (2 * pow(FertPar[1], 2)))
    
    if np.any(TempFunc(TemPar, i) >= FertPar[3]) and np.any(TempFunc(TemPar, i) <= FertPar[4]):
        FP = FP_1 * FP_2 * FP_3
    else:
       FP = 0
    return FP


# Definition of the Mortality rate function

def MortFunc(MortPar, TempPar, i):
    return MortPar[0] * (pow(TempFunc(TemPar, i) - MortPar[1], 2))  


# Definition of the random parameters


BriPar = BriPar_generator(a, T_L, T_M, m)
BioPar = BioPar_generator(a_1, h_b, P, var)
FertPar = FertPar_generator(gamma, sigma_T, T_star, T_low, T_max, var)
MortPar = MortPar_generator(a_2, T_opt, var)


# Definition of the time step and algorithm for PDE integration:
#    the variables time, age and biomass are iteration variables
#    the "delta" variables are the step increases

def dN(N, time, age, biomass, delta_i, delta_h, delta_s, Br, m, f, I_t, G):

    #SOL =  N[time, age, biomass] - Br * ((time * delta_i)/(age * delta_h)) * (N[time, age, biomass] - N[time, age - 1, biomass]) - G * ((time * delta_i)/(biomass * delta_s)) * (N[time, age, biomass] - N[time, age, biomass - 1]) + (time * delta_i) * m * N[time, age, biomass] + (time * delta_i) * f * I_t
    SOL =  N[time, age, biomass] - Br * (delta_i/delta_h) * (N[time, age, biomass] - N[time, age - 1, biomass]) - G * (delta_i/delta_s) * (N[time, age, biomass] - N[time, age, biomass - 1]) + delta_i * m * N[time, age, biomass] + delta_i * f * I_t
    return SOL


# Definition of the numerical scheme for PDE solution

def PDE_solution(i, h, s, delta_i, delta_h, delta_s, BriPar, BioPar, FertPar, MortPar, N_iter):

    # Open the storage files on the basis of the iteration

    name_file = 'Solutions_' + str(N_iter) + '.txt'
    res_txt = open(name_file, 'w')

    # Create the layers of the matrix (it is a cube!)
    
    # Firts layer, the plans (time)

    N_1 = np.zeros((i))[:,None, None] # time, age, biomass in order
    
    # Second layer, the rows (age)
    
    N_2 = np.zeros((h))[:, None]
    
    # Third layer, the columns (biomass)
    
    N_3 = np.zeros((s))
    
    # Set the arrays all together (create the matrix)
    
    N = N_1 + N_2 + N_3
    
    # Set the array for integration

    I_t_store = np.zeros((h))

    for time in range(i):
        for age in range(h):
            for biomass in range(s):

                if time == 0 and age == 0:

                    # Set initial and boundary conditions

                    N[0,:,:]= 10
                    N[:,0,:]= BriFunc(BriPar, TemPar, time*delta_i)
                    
                else:
                    m = MortFunc(MortPar, TemPar, time*delta_i)
                    f = FertFunc(FertPar, TemPar, time*delta_i)
                    Br = BriFunc(BriPar, TemPar, time*delta_i)
                    I_t_store[age] += N[time, age, biomass]
                    G = BioFunc(BioPar)

                    I_t = I_t_store[age]

                    N[time, age, biomass] = dN(N, time, age, biomass, delta_i, delta_h, delta_s, Br, m, f, I_t, G) + N[time, age, biomass]

                res_txt.write(str(time * delta_i) + '\t' + str(age * delta_h) + '\t' + str(biomass * delta_s) + '\t' + str(N[time, age, biomass]) + '\n')
    
    res_txt.close()













 

