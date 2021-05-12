"""===========================================================================


TITLE: 

AUTH:    Anthony JIMENEZ
DOB:     XX MMMM YYYY
LAST:    XX MMMM YYYY

DESCR:   Assignment 05 Re



=========================================================================="""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.special as sp
from scipy.optimize import shgo
from scipy.optimize import minimize

## EXAMPLE PROBLEM ASSUMPTIONS
SG = 0.75     ## DIM-LESS
p_pc = 667    ## PSIA
T_pc = 408    ## DEGREE RANKINE
p_ts = 2500   ## PSIA
T_ts = 35     ## DEGREE FAHRENHEIT
T_ws = 245    ## DEGREE FAHRENHEIT
L = 10000     ## FEET
n_cut = 2     ## NUMBER OF DISCRETIZATIONS ACROSS LENGTH
theta = 0     ## DEGREES
z_ts = 0.585  ## COMPUTED FROM DRANCHUK-ABOU-KASSEM CORRELATION
eps_stop = .1  ## STOPPING CRITERIA FOR CONVERGENCE

## OILFIELD UNIT CONSTANTS AND CONVERSIONS
R = 10.731557089016  ## GAS CONSTANT (FIELD UNITS)
M_air = 28.97        ## MOLAR MASS OF AIR (g/MOL)
grav_c = 144         ## GRAVITATIONAL CONVERSION CONSTANT (FIELD UNITS)
alpha_T = 459.67      ## TEMPERATURE CONVERSION CONSTANT (F + alpha_T ---> R)
alpha_c = M_air / R / grav_c


def rhs_const_func(SG, L, theta):
    return alpha_c * SG * L * np.cos(theta)
def lhs_integrand_func(z, T, p):
    T = T + alpha_T
    return T * z / p
def DAK_z_factor_func(z, T_mid, p_mid):
    a = np.array([0.3265, -1.0700, -0.5339, 0.01569, -0.05165,
                  0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210])
    Tr = (T_mid + alpha_T) / T_pc
    pr = p_mid / p_pc
    dens_red = 0.27 * pr / z / Tr
    term1 = (a[0] + a[1]/Tr + a[2]/Tr**3 + a[3]/Tr**4 + a[4]/Tr**5) * dens_red
    term2 = (a[5] + a[6]/Tr + a[7]/Tr**2) * dens_red**2
    term3 = (a[6]/Tr + a[7]/Tr**2) * a[8] * dens_red**5
    term4 = a[9] * (1 + a[10] * dens_red**2) * (dens_red**2 / Tr**3) * np.exp(-a[10] * dens_red**2)
    obj_func = (1 + term1 + term2 - term3 + term4 - z)**2
    return obj_func
def iterate_segment_func(eps_stop, SG, L, theta, n_cut, z, T_top, p_top, i_seg):
    count = 0
    rel_error = 100
    const_rhs = rhs_const_func(SG, L, theta)
    T_mid = T_top + (T_ws - T_ts)/L * (L/(n_cut+1)*i_seg)
    while (rel_error > eps_stop) == True:
        if count == 0:
            I_top = lhs_integrand_func(z, T_top, p_top)
            I_mid = I_top
            p_mid = p_top + const_rhs / (I_mid + I_top)
            count += 1
        else:
            p_mid_old = p_mid
            z_mid = minimize(DAK_z_factor_func, z, args= (T_mid, p_mid))
            z_mid = z_mid.x
            I_mid = lhs_integrand_func(z_mid, T_mid, p_mid)
            p_mid = p_top + const_rhs / (I_mid + I_top)
            rel_error = np.abs((p_mid - p_mid_old) / p_mid_old) * 100
            count += 1
        if count > 100:
            print('FAILED TO CONVERGE')
            exit
    return p_mid, z_mid, T_mid, rel_error
'''
def base_func(eps_stop, SG, L, theta, n_cut, z, T_top, p_top):
    p_arr = np.zeros([n_cut+1])
    for i in range(n_cut+1):
        p_arr[i], z, T_top, rel_error = iterate_segment_func(eps_stop, SG, L, theta, n_cut, z, T_top, p_top)
        p_top = p_arr[i]
    return p_arr

p_grad = base_func(eps_stop, SG, L, theta, n_cut, z_ts, T_ts, p_ts)
'''

p_arr = np.zeros([n_cut+2])
for i in range(n_cut+2):
    if i == 0:
        p_arr[i] = p_ts
    else:
        p_arr[i], z, T_top, rel_error = iterate_segment_func(eps_stop, SG, L, theta, n_cut, z_ts, T_ts, p_ts, i)
        T_ts = T_top
        z_ts = z
        p_ts = p_arr[i]
        print(z_ts)

#%%

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.special as sp
from scipy.optimize import shgo
from scipy.optimize import minimize

## EXAMPLE PROBLEM ASSUMPTIONS
SG = 0.75       ## DIM-LESS
p_pc = 667      ## PSIA
T_pc = 408      ## DEGREE RANKINE
p_ts = 2500     ## PSIA
T_ts = 35       ## DEGREE FAHRENHEIT
T_ws = 245      ## DEGREE FAHRENHEIT
L = 10000       ## FEET
n_cut = 1       ## NUMBER OF DISCRETIZATIONS ACROSS LENGTH
theta = 0       ## DEGREES
z_ts = 0.585    ## COMPUTED FROM DRANCHUK-ABOU-KASSEM CORRELATION
eps_stop = 0.1  ## STOPPING CRITERIA FOR CONVERGENCE

## OILFIELD UNIT CONSTANTS AND CONVERSIONS
R = 10.731557089016  ## GAS CONSTANT (FIELD UNITS)
M_air = 28.97        ## MOLAR MASS OF AIR (g/MOL)
grav_c = 144         ## GRAVITATIONAL CONVERSION CONSTANT (FIELD UNITS)
alpha_T = 459.67      ## TEMPERATURE CONVERSION CONSTANT (F + alpha_T ---> R)
alpha_c = M_air / R / grav_c


def rhs_const_func(SG, L, theta):
    return alpha_c * SG * L * np.cos(theta)
def lhs_integrand_func(z, T, p):
    T = T + alpha_T
    return T * z / p
def DAK_z_factor_func(z, T_mid, p_mid):
    a = np.array([0.3265, -1.0700, -0.5339, 0.01569, -0.05165,
                  0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210])
    Tr = (T_mid + alpha_T) / T_pc
    pr = p_mid / p_pc
    dens_red = 0.27 * pr / z / Tr
    term1 = (a[0] + a[1]/Tr + a[2]/(Tr**3) + a[3]/(Tr**4) + a[4]/(Tr**5)) * dens_red
    term2 = (a[5] + a[6]/Tr + a[7]/(Tr**2)) * dens_red**2
    term3 = (a[6]/Tr + a[7]/(Tr**2)) * a[8] * dens_red**5
    term4 = a[9] * (1 + a[10] * (dens_red**2)) * ((dens_red**2) / (Tr**3)) * np.exp(-a[10] * (dens_red**2))
    obj_func = (1 + term1 + term2 - term3 + term4 - z)**2
    return obj_func
def iterate_segment_func(eps_stop, SG, L, theta, n_cut, z, T_top, p_top, T_0, i_seg):
    count = 0
    rel_error = 100
    const_rhs = rhs_const_func(SG, L, theta) / n_cut
    T_mid = T_top + (T_ws - T_top)/L * (L/(n_cut+1)*i_seg)
    I_top = lhs_integrand_func(z, T_top, p_top)
    while (rel_error > eps_stop) == True:
        if count == 0:
            I_mid = I_top
            p_mid = p_top + const_rhs / (I_top + I_mid)
            count += 1
        else:
            p_mid_old = p_mid
            z = minimize(DAK_z_factor_func, z, args= (T_mid, p_mid))
            z = z.x
            I_mid = lhs_integrand_func(z, T_mid, p_mid)
            p_mid = p_top + const_rhs / (I_top + I_mid)
            rel_error = np.abs((p_mid - p_mid_old) / p_mid_old) * 100
            count += 1
    if count > 100:
        print('FAILED TO CONVERGE')
        exit
    return p_mid, z


p_arr = np.zeros([n_cut+3])
T_0 = T_ts
for i in range(n_cut + 3):
    print(p_ts, z_ts, T_ts)
    if i == 0:
        p_arr[i] = p_ts
    else:
        p_arr[i], z = iterate_segment_func(eps_stop, SG, L, theta, n_cut, z_ts, T_ts, p_ts, T_0, i)
        z_ts = z
        p_ts = p_arr[i]
        








dL = np.linspace(0, n_cut+1, n_cut+2) * (L / (n_cut+1))
#z = minimize(DAK_z_factor_func, 0.5, args= (140, 3086))
#print(z.x*(140+459.67) / 3086)
#print(z.x)




#%%
plt_lw = 2
fig, ax = plt.subplots(figsize= [7,7], dpi= 150)
ax.plot(x, y, c= 'k', lw= plt_lw, label= 'NA')
ax.set(xlabel= 'NA', ylabel= 'NA', xlim= (0,1), ylim= (0,1))
ax.minorticks_on()
plt.title('NA')
plt.legend(framealpha= 1)
plt.show()

