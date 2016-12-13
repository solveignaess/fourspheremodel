import os
import numpy as np
from parameters_wm import *

def conductivity(sigma_skull):
    s12 = sigma_brain / sigma_csf
    s23 = sigma_csf / sigma_skull
    s34 = sigma_skull / sigma_scalp
    return s12, s23, s34

# scalp_rad = scalp_rad - rad_tol

rz = dipole_loc
rz1 = rz / brain_rad
r12 = brain_rad / csftop_rad
r23 = csftop_rad / skull_rad
r34 = skull_rad / scalp_rad

r1z = 1. / rz1
r21 = 1. / r12
r32 = 1. / r23
r43 = 1. / r34

def V(n):
    k = (n+1.) / n
    Factor = ( ( r34**n - (r43**(n+1)) ) / ( (k*(r34**n)) + (r43**(n+1)) ) )
    num = (s34*k) - Factor
    den = s34 + Factor
    return (num / den)

def Y(n):
    k = n / (n+1.)
    Factor = ( (r23**n) * ( (k - V(n)*(r32**(n+1))) / (r23**n + V(n)*(r32**(n+1)) )))
    num = (s23*k) - Factor
    den = s23 + Factor
    return (num / den)

def Z(n):
    k = (n+1.) / n
    num = (r12**n - k*Y(n)*(r21**(n+1)) ) / (r12**n + Y(n)*(r21**(n+1)))
    return num

def A1(n):
    num = (rz1**(n-1))*((n*Z(n)) + s12*(n+1.))
    den = (s12 - Z(n))*4*np.pi*(brain_rad**2)
    return num / den

def A2(n):
    num = A1(n) + (n*(rz1**(n-1)) / (4*np.pi*(brain_rad**2)))
    den = (Y(n)*(r21**(n+1))) + r12**n
    return num / den

def B2(n):
    return A2(n)*Y(n)

def A3(n):
    num = A2(n) + B2(n)
    den = r23**n + (V(n)*(r32**(n+1)))
    return num / den

def B3(n):
    return A3(n)*V(n)

def A4(n):
    k = (n+1.) / n
    num = A3(n) + B3(n)
    den = (k*(r34**n)) + (r43**(n+1))
    return k*(num / den)

def B4(n):
    return A4(n)* (n / (n+1.))

def H(n):
    return (A4(n) + B4(n))

I = 1.
d = .1
n = np.arange(1, 10000)

s12, s23, s34 = conductivity(sigma_skull20)
coef = H(n)
coef = np.insert(coef, 0, 0)
Lprod = np.polynomial.legendre.Legendre(coef)
phi_20 = I*d*Lprod(np.cos(np.deg2rad(theta)))  / sigma_brain
np.save(os.path.join('results', 'phi_20_98.npy'), phi_20)

s12, s23, s34 = conductivity(sigma_skull40)
coef = H(n)
coef = np.insert(coef, 0, 0)
Lprod = np.polynomial.legendre.Legendre(coef)
phi_40 = I*d*Lprod(np.cos(np.deg2rad(theta))) / sigma_brain
np.save(os.path.join('results', 'phi_40_98.npy'), phi_40)

s12, s23, s34 = conductivity(sigma_skull80)
coef = H(n)
coef = np.insert(coef, 0, 0)
Lprod = np.polynomial.legendre.Legendre(coef)
phi_80 = I*d*Lprod(np.cos(np.deg2rad(theta))) / sigma_brain
np.save(os.path.join('results', 'phi_80_98.npy'), phi_80)
