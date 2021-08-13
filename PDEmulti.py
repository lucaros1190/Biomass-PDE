
# Python script to manage the multiproces of PDE-biomass.py script
# Created by Luca Rossini
# e-mail: luca.rossini@unitus.it
# Last update 17 July 2021


# List of import

import matplotlib.pyplot as plt
import matplotlib
import multiprocessing as mp
import os
import time
from PDEbiomass import *
from Parameters import *
from matplotlib.pyplot import cm


# Function which generates random parameters and PDE solution

def PDE_random(a, T_L, T_M, m, a_1, h_b, P, var, gamma, sigma_T, T_star, T_low, T_max, a_2, T_opt, N_iter, i, h, s, delta_i, delta_h, delta_s):
    BriPar = BriPar_generator(a, T_L, T_M, m)
    BioPar = BioPar_generator(a_1, h_b, P, var)
    FertPar = FertPar_generator(gamma, sigma_T, T_star, T_low, T_max, var)
    MortPar = MortPar_generator(a_2, T_opt, var)
    PDE_solution(i, h, s, delta_i, delta_h, delta_s, BriPar, BioPar, FertPar, MortPar, N_iter)


# Manage the multiprocess for equation solving

mp.Pool(processes = N_iter)

for processes in range(N_iter):
    MultiProcess_PDE = mp.Process(target=PDE_random, args = (a, T_L, T_M, m, a_1, h_b, P, var, gamma, sigma_T, T_star, T_low, T_max, a_2, T_opt, processes, i, h, s, delta_i, delta_h, delta_s))

    MultiProcess_PDE.start()

time.sleep(sleep)


# Functions which make the final plots

# Biomass vs population

def bio_vs_pop(N_iter):

    plt.figure(1)

    color=iter(cm.rainbow(np.linspace(0,1,N_iter-1)))

    for plot in range(N_iter-1):

        name_file = 'Solutions_' + str(plot+1) + '.txt'
        data = pd.read_csv(name_file, sep="\t", header=0, engine='python')

        data.columns = ["Time", "Age", "Biomass", "Pop_density"]
        bio = data['Biomass']
        pop = data['Pop_density']
    
        c = next(color)
        plt.plot(bio, pop, '.', label = 'Population', c = c, markersize = 0.25, rasterized=True)
        plt.xlabel('Biomass (Biomass classes)')
        plt.ylabel('Population density (N of individuals)')
        plt.axis('tight')


# Age vs population

def age_vs_pop(N_iter):

    plt.figure(2)

    color=iter(cm.rainbow(np.linspace(0,1,N_iter-1)))

    for plot in range(N_iter-1):

        name_file = 'Solutions_' + str(plot+1) + '.txt'
        data = pd.read_csv(name_file, sep="\t", header=0, engine='python')

        data.columns = ["Time", "Age", "Biomass", "Pop_density"]
        ag = data['Age']
        pop = data['Pop_density']
    
        c = next(color)
        plt.plot(ag, pop, '.', label = 'Population', c = c, markersize = 0.25, rasterized=True)
        plt.xlabel('Age (Stage number)')
        plt.ylabel('Population density (N of individuals)')
        plt.axis('tight')


# Time vs population

def time_vs_pop(N_iter):

    plt.figure(3)

    color=iter(cm.rainbow(np.linspace(0,1,N_iter-1)))

    for plot in range(N_iter-1):

        name_file = 'Solutions_' + str(plot+1) + '.txt'
        data = pd.read_csv(name_file, sep="\t", header=0, engine='python')

        data.columns = ["Time", "Age", "Biomass", "Pop_density"]
        t = data['Time']
        pop = data['Pop_density']
    
        c = next(color)
        plt.plot(t, pop, '.', label = 'Population', c = c, markersize = 0.25, rasterized=True)
        plt.xlabel('Time (days)')
        plt.ylabel('Population density (N of individuals)')
        plt.axis('tight')


# Time vs age

def time_vs_age(N_iter):

    plt.figure(4)

    color=iter(cm.rainbow(np.linspace(0,1,N_iter-1)))

    for plot in range(N_iter-1):

        name_file = 'Solutions_' + str(plot+1) + '.txt'
        data = pd.read_csv(name_file, sep="\t", header=0, engine='python')

        data.columns = ["Time", "Age", "Biomass", "Pop_density"]
        ag = data['Age']
        time = data['Time']
    
        c = next(color)
        plt.plot(time, ag, '.', label = 'Population', c = c, markersize = 0.25, rasterized=True)
        plt.xlabel('Time (days)')
        plt.ylabel('Age (Stage number)')
        plt.axis('tight')


# Time vs biomass

def time_vs_biomass(N_iter):

    plt.figure(5)

    color=iter(cm.rainbow(np.linspace(0,1,N_iter-1)))

    for plot in range(N_iter-1):

        name_file = 'Solutions_' + str(plot+1) + '.txt'
        data = pd.read_csv(name_file, sep="\t", header=0, engine='python')

        data.columns = ["Time", "Age", "Biomass", "Pop_density"]
        biomass = data['Biomass']
        time = data['Time']
    
        c = next(color)
        plt.plot(time, biomass, '.', label = 'Population', c = c, markersize = 0.25, rasterized=True)
        plt.xlabel('Time (days)')
        plt.ylabel('Biomass (Biomass classes)')
        plt.axis('tight')


# Make the 3D-plot


def plot3d(N_iter):

    plt.figure(6)

    for plot in range(N_iter-1):

        name_file = 'Solutions_' + str(plot+1) + '.txt'
        data = pd.read_csv(name_file, sep="\t", header=0, engine='python')

        data.columns = ["Time", "Age", "Biomass", "Pop_density"]
        t = data['Time']
        ag = data['Age']
        bio = data['Biomass']
        pop = data['Pop_density']
    
        ax = plt.axes(projection='3d')
        ax.plot3D(t, bio, pop)
        plt.axis('tight')


# Plot the results using the functions defined below

Biomass_plot = bio_vs_pop(N_iter)
Age_plot = age_vs_pop(N_iter)
Time_plot = time_vs_pop(N_iter)
Time_age_plot = time_vs_age(N_iter)
Time_biomass_plot = time_vs_biomass(N_iter)

TreDPlot = plot3d(N_iter)

plt.show()



