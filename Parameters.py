import multiprocessing as mp

# Define variables

N_iter = mp.cpu_count() 
sleep = 30

# Random variation of the parameters and seeding

var = 0.05

# Parameter values and step for plotting:

i = 3600
h = 8
s = 10

delta_i = 0.1 # Time step increase
delta_h = 1 # Age step increase
delta_s = 1 # Biomass step increase

# Briere 

a = 3.17 * pow(10, -5)
T_L = 9.00
T_M = 38.99
m = 3.07 

# Biomass

a_1 = 0.00001
h_b = (1/240)
P = 200

# Fertility

gamma = 100.0
sigma_T = 4.0
T_star = 33.61
T_low = 13.0
T_max = 40.0

# Mortality

a_2 = 2 * pow(10, -5)
T_opt = 35.66

# Temperature - Not subject on random variations!

T_0 = 30.0
A = 1.0
Period = 120 # Unit in days
TemPar = [T_0, A, Period]


