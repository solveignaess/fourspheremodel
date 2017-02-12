import os
import numpy as np
from scipy.special import lpmv
import parameters as params


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
    den = (s12 - Z(n))*4*np.pi*(params.brain_rad**2)
    return num / den


def A2(n):
    num = A1(n) + (n*(rz1**(n-1)) / (4*np.pi*(params.brain_rad**2)))
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


def H(n, r_ele=params.scalp_rad):
    if r_ele < params.brain_rad:
        T1 = ((r_ele / params.brain_rad)**n) * A1(n)
        T2 = ((rz / r_ele)**(n + 1))
    elif r_ele < params.csftop_rad:
        T1 = ((r_ele / params.csftop_rad)**n) * A2(n)
        T2 = ((params.csftop_rad / r_ele)**(n + 1)) * B2(n)
    elif r_ele < params.skull_rad:
        T1 = ((r_ele / params.skull_rad)**n) * A3(n)
        T2 = ((params.skull_rad / r_ele)**(n + 1)) * B3(n)
    elif r_ele <= params.scalp_rad:
        T1 = ((r_ele / params.scalp_rad)**n) * A4(n)
        T2 = ((params.scalp_rad / r_ele)**(n + 1)) * B4(n)
    else:
        print "Invalid electrode position"
        return
    return (T1 + T2)


def decompose_dipole(I):
    P = np.array([np.array(src_pos) * I - np.array(snk_pos) * I])
    dp_loc = (np.array(src_pos) + np.array(snk_pos)) / 2.
    dist_dp = np.linalg.norm(dp_loc)
    dp_rad = (np.dot(P, dp_loc) / dist_dp) * (dp_loc / dist_dp)
    dp_tan = P - dp_rad
    sign_rad = np.sign(np.dot(P, dp_loc))
    dp_rad = sign_rad * np.linalg.norm(dp_rad)
    dp_tan = np.linalg.norm(dp_tan)
    return dp_rad, dp_tan


def conductivity(sigma_skull):
    s12 = params.sigma_brain / params.sigma_csf
    s23 = params.sigma_csf / sigma_skull
    s34 = sigma_skull / params.sigma_scalp
    return s12, s23, s34


def compute_phi(s12, s23, s34, dp_rad, dp_tan):
    coef = H(n)
    cos_theta = np.cos(params.theta_r)

    # radial
    n_coef = coef
    rad_coef = np.insert(n_coef, 0, 0)
    Lprod = np.polynomial.legendre.Legendre(rad_coef)
    Lfactor_rad = Lprod(cos_theta)
    rad_phi = dp_rad * Lfactor_rad

    # # #tangential
    # Lfuncprod = []
    # for tt in range(params.theta_r.size):
    #     Lfuncprod.append(np.sum([C * lpmv(1, P, cos_theta[tt])
    #                              for C, P in zip(coef, n)]))
    # tan_phi = -1 * dp_tan * np.sin(params.phi_angle_r) * np.array(Lfuncprod)
    return (rad_phi ) / (params.sigma_brain)


# scalp_rad = scalp_rad - rad_tol
rz = params.dipole_loc
rz1 = rz / params.brain_rad
r12 = params.brain_rad / params.csftop_rad
r23 = params.csftop_rad / params.skull_rad
r34 = params.skull_rad / params.scalp_rad

r1z = 1. / rz1
r21 = 1. / r12
r32 = 1. / r23
r43 = 1. / r34

I = 1.
n = np.arange(1, 100)

dipole = params.dipole_list[0]  # 'rad_dipole'

print 'WARNING: These results are for comparision only!'
print 'Please use the correct formulation instead'
print 'Now computing for dipole using Srinivasan98: ', dipole['name']
src_pos = dipole['src_pos']
snk_pos = dipole['snk_pos']
dp_rad, dp_tan = decompose_dipole(I)
print 'Not evaluating the tangential component'

s12, s23, s34 = conductivity(params.sigma_skull20)
phi_20 = compute_phi(s12, s23, s34, dp_rad, dp_tan)

s12, s23, s34 = conductivity(params.sigma_skull40)
phi_40 = compute_phi(s12, s23, s34, dp_rad, dp_tan)

s12, s23, s34 = conductivity(params.sigma_skull80)
phi_80 = compute_phi(s12, s23, s34, dp_rad, dp_tan)

s12 = s23 = s34 = 1.
phi_lim = compute_phi(s12, s23, s34, dp_rad, dp_tan)

f = open(os.path.join('results',
                      'Analytical_Sri98_' + dipole['name'] + '.npz'), 'w')
np.savez(f, phi_20=phi_20, phi_40=phi_40, phi_80=phi_80, phi_lim=phi_lim)
f.close()