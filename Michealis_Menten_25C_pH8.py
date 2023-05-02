#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 17:47:55 2023

@author: gradycorkum
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 17:14:23 2023

@author: gradycorkum
"""

# Imports
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#%% Simulate MM kinetics of mutant enzyme at 25c ph 8
# Constants
K_cat = 136.21
K_M = 19.28
# Variables
E_0 = 10e-3
# ODE
def MM_dynamics(t,y,K_cat,K_M,E_0):
    # y[0] = B product
    # y[1] = A reactant
    dydt = np.zeros(2)
    v = K_cat*E_0*(y[1]/(K_M+y[1]))
    dydt[0] = v
    dydt[1] = -v
    return dydt


# Initial Conditions
y0 = [0,100]
# time span
t = np.linspace(0,50,100)
tspan = [t[0],t[-1]]

ode_sol = solve_ivp(lambda t,y:MM_dynamics(t,y,K_cat,K_M,E_0),tspan,y0,t_eval=t)

plt.figure(1)
plt.title("Simulate MM Dynamics of Recombinant @25C pH 8")
plt.plot(t,ode_sol.y[0],'k-')
plt.plot(t,ode_sol.y[1],'r--')
plt.xlabel('time [s]')
plt.ylabel('concentration [mM]')
plt.legend(['product','reactant']);
plt.tight_layout()

#%% Simulate initial reactant concentration
plt.figure(2)

r_conc = 100

# time span
t = np.linspace(0,80,10000)
tspan = [t[0],t[-1]]

# Initial Conditions
y0 = [0,r_conc]
ode_sol = solve_ivp(lambda t,y:MM_dynamics(t,y,K_cat,K_M,E_0),tspan,y0,t_eval=t)
r = ode_sol.y[1] #reactant concentration

# Plot to make sure simulation has run for long enough
plt.title("Simulate MM Dynamics (half lives)")
plt.plot(t,ode_sol.y[0],'k-')
plt.plot(t,ode_sol.y[1],'r--')
plt.xlabel('time [s]')
plt.ylabel('concentration [mM]')
plt.legend(['product','reactant']);
plt.tight_layout()

# Calculate time from 100% to 50%
t100_50 = t[np.argmin(np.abs(r-r[0]*.5))]
print(str(r_conc),'mM initial reactant, time from 100% to 50%: ',np.round(t100_50,2),'seconds')

# Calculate time from 50% to 25%
t50_25 = t[np.argmin(np.abs(r-r[0]*0.25))]-t[np.argmin(np.abs(r-r[0]*.5))]
print(str(r_conc),'mM initial reactant, time from 50% to 25%:',np.round(t50_25,2),'seconds')


#%%
plt.figure(3)
# Same parameters as before
K_cat = K_cat
K_M = K_M
# Variables
E_0 = 10e-3
# Range of reactant concentrations
A = np.linspace(0,80,1000)
# Calculate rate from MM equation
v = K_cat*E_0*(A/(K_M+A))
V_Max = E_0 * K_cat


# plot
plt.plot(A,v)
#plot Vmax
plt.plot([K_M,K_M ],[0,1.75])
plt.plot([0,80], [E_0 * K_cat, E_0 * K_cat])
plt.xlabel('[A] [mM]')
plt.ylabel('V [mM/s]')
plt.title("Reaction Rate")
plt.tight_layout()
plt.legend(["Rate","K_M","VMax"]);



#%% Lineweaver Burke Plot
plt.figure(4)
plt.title('Lineweaver Burke Plot')
# plot
plt.plot(1/A[1:],1/v[1:]) #we won't plot the first point to avoid a divide by 0 warning
plt.xlabel('1/[A]')
plt.ylabel('1/V')



