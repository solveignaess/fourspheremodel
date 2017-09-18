# -*- coding: utf-8 -*-
from __future__ import division
import os
from scipy.special import lpmv
import numpy as np
from plotting_convention import simplify_axes, mark_subplots
import sys
reload(sys)
sys.setdefaultencoding('UTF-8')
class CalcPotential4Sphere:
    """ Import all data needed for calculating the extracellular potential Phi at location r
    from a current dipole moment p at location rz. Also decompose p into a radial and a tangential part.
    """

    def __init__(self, radii, sigmas, r, rz, source=None):
        # self.k1 = 1E3  # from mV to muV
        self.k1 = 1E6  # from mV to nV
        self.r1 = radii[0]
        self.r2 = radii[1]
        self.r3 = radii[2]
        self.r4 = radii[3]
        # print self.r1, self.r2, self.r3, self.r4

        self.sigma1 = sigmas[0]
        self.sigma2 = sigmas[1]
        self.sigma3 = sigmas[2]
        self.sigma4 = sigmas[3]

        self.r12 = self.r1 / self.r2
        self.r21 = self.r2 / self.r1
        self.r23 = self.r2 / self.r3
        self.r32 = self.r3 / self.r2
        self.r34 = self.r3 / self.r4
        self.r43 = self.r4 / self.r3

        self.sigma12 = self.sigma1 / self.sigma2
        self.sigma21 = self.sigma2 / self.sigma1
        self.sigma23 = self.sigma2 / self.sigma3
        self.sigma32 = self.sigma3 / self.sigma2
        self.sigma34 = self.sigma3 / self.sigma4
        self.sigma43 = self.sigma4 / self.sigma3

        self.rxyz = r
        self.rzloc = rz
        self.rz = np.sqrt(np.sum(rz ** 2))
        self.rz1 = self.rz / self.r1
        self.r = np.sqrt(np.sum(r ** 2))

        self.source = source

    def calc_potential(self, p):
        """
        Parameters
        __________
        p : ndarray [1E-15 Am]
            Array containing the current dipole moment for all timesteps in
            the x-, y- and z-direction.

        Returns
        _______
        potential : ndarray [nV]
            Array containing the electric potential at point self.r.

        """
        # print 'calc potential'
        p_rad, p_tan = self.decompose_dipole(p)

        pot_rad = self.calc_rad_potential(p_rad)
        pot_tan = self.calc_tan_potential(p_tan)

        pot_tot = pot_rad + pot_tan
        # print p_rad, p_tan
        return pot_tot

    def decompose_dipole(self, p):
        """Decompose current dipole moment vector in radial and tangential terms
        Parameters
        __________
        p : ndarray [1E-15 Am]
            Array containing the current dipole moment for all timesteps in
            the x-, y- and z-direction.

        Returns:
        ________
        p_rad : ndarray [1E-15 Am]
                Radial part of p, parallel to self.rz
        p_tan : ndarray [1E-15 Am]
                Tangential part of p, orthogonal to self.rz
        """
        p_rad = np.array([np.dot(p, self.rzloc)/self.rz ** 2 * self.rzloc])
        p_tan = p - p_rad

        return p_rad, p_tan

    def calc_rad_potential(self, p_rad):
        """Return potential from radial dipole p_rad at location rz measured at r
        Parameters
        __________
        P : ndarray [1E-15 Am]
            Array containing the current dipole moment for all timesteps in
            the x-, y- and z-direction.

        Returns
        _______
        potential : ndarray [nV]
            Array containing the current dipole moment at point r.
        """

        # print 'calc rad pot'
        p_tot = np.linalg.norm(p_rad)
        theta = self.calc_theta()
        s_vector = self.sign_rad_dipole(p_rad)
        s_vector = s_vector.reshape(len(p_rad), 1)
        phi_const = s_vector * p_tot / (4 * np.pi * self.sigma1 * self.rz ** 2) * self.k1
        if self.r <= self.r1:
            # print 'rad brain', self.r, self.r1
            potential_nterms = self.potential_brain_rad(theta)
        elif self.r <= self.r2:
            # print 'rad csf'
            potential_nterms = self.potential_csf_rad(theta)
        elif self.r <= self.r3:
            # print 'rad skull'
            potential_nterms = self.potential_skull_rad(theta)
        elif self.r <= (self.r4):
            # print 'rad scalp'
            potential_nterms = self.potential_scalp_rad(theta)
        elif self.r <= (self.r4+1E-6):
            self.r = self.r4
            potential_nterms = self.potential_scalp_rad(theta)
        else:
            potential_nterms = np.nan
            raise ValueError('Electrode located outside head model. Maximum r = %s µm.', self.r4, '\n your r = ', self.r)
        potential = phi_const * potential_nterms
        return potential

    def calc_tan_potential(self, p_tan):
        """Return potential from tangential dipole P at location rz measured at r
        Parameters
        __________
        p : ndarray [1E-15 Am]
            Array containing the current dipole moment for all timesteps in
            the x-, y- and z-direction.

        Returns
        _______
        potential : ndarray [nV]
            Array containing the current dipole moment at point r.
        """

        # print 'calc tan pot'
        theta = self.calc_theta()
        phi = self.calc_phi(p_tan)
        p_tot = np.linalg.norm(p_tan, axis=1)
        phi_hom = - p_tot / (4 * np.pi * self.sigma1 * self.rz ** 2) * np.sin(phi) * self.k1
        if self.r <= self.r1:
            # print 'tan brain'
            potential_nterms = self.potential_brain_tan(theta)
        elif self.r <= self.r2:
            # print 'tan csf'
            potential_nterms = self.potential_csf_tan(theta)
        elif self.r <= self.r3:
            # print 'tan skull'
            potential_nterms = self.potential_skull_tan(theta)
        elif self.r <= self.r4:
            # print 'tan scalp'
            potential_nterms = self.potential_scalp_tan(theta)
        elif self.r <= (self.r4+1E-6):
            self.r = self.r4
            potential_nterms = self.potential_scalp_rad(theta)
        else:
            potential_nterms = np.nan
            raise ValueError('Electrode located outside head model. Maximum r = %s µm.', self.r4)

        potential = phi_hom * potential_nterms
        return potential

    def calc_theta(self):
        """Calculate theta and phi for radial dipole"""
        cos_theta = np.dot(self.rzloc, self.rxyz) / (np.linalg.norm(self.rxyz) * np.linalg.norm(self.rzloc))
        cos_theta = np.nan_to_num(cos_theta)
        theta = np.arccos(cos_theta)
        return theta

    def calc_phi(self, p):
        """Calculate phi-angles for tangential dipole."""
        proj_rxyz_rz = np.dot(self.rxyz, self.rzloc) / np.sum(self.rzloc ** 2) * self.rzloc
        rxy = self.rxyz - proj_rxyz_rz
        rx = np.cross(p, self.rzloc)/(np.linalg.norm(p) * self.rz)
        cos_phi = np.dot(rx, rxy) / (np.linalg.norm(rxy) * np.linalg.norm(rx))
        cos_phi = np.nan_to_num(cos_phi)
        phi_wrong_sign = np.arccos(cos_phi)

        phi = phi_wrong_sign
        for i in range(len(p)):
            sign_test = np.dot(p, self.rxyz)/(self.r*np.linalg.norm(p))
            if sign_test < 0:
                phi[i] = 2*np.pi - phi_wrong_sign[i]
        return phi


    def sign_rad_dipole(self, p):
        """Flip radial dipole pointing inwards, and add a -1
        to the s-vector, so that the potential can be
        calculated as if the dipole were pointing outwards,
        and later be multiplied by -1 to get the right potential."""
        sign_vector = np.ones((len(p), 1))
        for i in range(len(p)):
            radial_test = np.dot(p, self.rzloc) / (np.linalg.norm(p) * self.rz)
            if np.abs(radial_test + 1) < 10 ** -8:
                # print p[i], np.dot(p, self.rzloc)
                sign_vector[i] = -1.
        return sign_vector

    def potential_brain_rad(self, theta):
        """Calculate sum with constants and legendres
        Parameters
        __________
        theta : ndarray [radians]
                Array of angles between electrode location and
                dipole location vectors.
        Returns
        _______
        pot_sum : float
                Sum factor containing brain constants, dipole and
                electrode locations and legendre polynomials.
        """
        n = np.arange(1, 100)
        c1n = self.calc_c1n(n)
        consts = n*(c1n * (self.r / self.r1) ** n + (self.rz / self.r) ** (n + 1))
        consts = np.insert(consts, 0, 0)
        leg_consts = np.polynomial.legendre.Legendre(consts)
        pot_sum = leg_consts(np.cos(theta))
        return pot_sum

    def potential_csf_rad(self, theta):
        n = np.arange(1,100)
        c2n = self.calc_c2n(n)
        d2n = self.calc_d2n(n, c2n)
        consts = n*(c2n * (self.r / self.r2) ** n + d2n * (self.r2 / self.r) ** (n + 1))
        consts = np.insert(consts, 0, 0)
        leg_consts = np.polynomial.legendre.Legendre(consts)
        pot_sum = leg_consts(np.cos(theta))
        return pot_sum

    def potential_skull_rad(self, theta):
        n = np.arange(1,100)
        c3n = self.calc_c3n(n)
        d3n = self.calc_d3n(n, c3n)
        consts = n*(c3n * (self.r / self.r3) ** n + d3n * (self.r3 / self.r) ** (n + 1))
        consts = np.insert(consts, 0, 0)
        leg_consts = np.polynomial.legendre.Legendre(consts)
        pot_sum = leg_consts(np.cos(theta))
        return pot_sum

    def potential_scalp_rad(self, theta):
        """Calculate potential in scalp from radial dipole.
        """
        n = np.arange(1,100)
        c4n = self.calc_c4n(n)
        d4n = self.calc_d4n(n, c4n)
        consts = n*(c4n * (self.r / self.r4) ** n + d4n * (self.r4 / self.r) ** (n + 1))
        consts = np.insert(consts, 0, 0)
        leg_consts = np.polynomial.legendre.Legendre(consts)
        pot_sum = leg_consts(np.cos(theta))
        return pot_sum

    def potential_brain_tan(self, theta):
        """Calculate sum with constants and legendres
        Parameters
        __________
        theta : ndarray [radians]
                Array of angles between electrode location and
                dipole location vectors.
        p :     ndarray [1Ee-15 Am]
                Array containing the current dipole moment for all timesteps in
                the x-, y- and z-direction.
        Returns
        _______
        pot_sum : float
                Sum factor containing brain constants, dipole and
                electrode locations and legendre polynomials.
        """
        n = np.arange(1,100)
        c1n = self.calc_c1n(n)
        consts = (c1n * (self.r / self.r1) ** n + (self.rz / self.r) ** (n + 1))
        pot_sum = np.sum([c*lpmv(1, i, np.cos(theta)) for c,i in zip(consts,n)])
        return pot_sum

    def potential_csf_tan(self, theta):
        """Calculate sum with constants and legendres
        Parameters
        __________
        theta : ndarray [radians]
                Array of angles between electrode location and
                dipole location vectors.
        p :     ndarray [1Ee-15 Am]
                Array containing the current dipole moment for all timesteps in
                the x-, y- and z-direction.
        Returns
        _______
        pot_sum : float
                Sum factor containing brain constants, dipole and
                electrode locations and legendre polynomials.
        """
        n = np.arange(1,100)
        c2n = self.calc_c2n(n)
        d2n = self.calc_d2n(n, c2n)
        consts = c2n*(self.r/self.r2)**n + d2n*(self.r2/self.r)**(n+1)
        pot_sum = np.sum([c*lpmv(1, i, np.cos(theta)) for c,i in zip(consts,n)])
        return pot_sum

    def potential_skull_tan(self, theta):
        """Calculate sum with constants and legendres
        Parameters
        __________
        theta : ndarray [radians]
                Array of angles between electrode location and
                dipole location vectors.
        p :     ndarray [1Ee-15 Am]
                Array containing the current dipole moment for all timesteps in
                the x-, y- and z-direction.
        Returns
        _______
        pot_sum : float
                Sum factor containing brain constants, dipole and
                electrode locations and legendre polynomials.
        """
        n = np.arange(1,100)
        c3n = self.calc_c3n(n)
        d3n = self.calc_d3n(n, c3n)
        consts = c3n*(self.r/self.r3)**n + d3n*(self.r3/self.r)**(n+1)
        pot_sum = np.sum([c*lpmv(1, i, np.cos(theta)) for c,i in zip(consts,n)])
        return pot_sum

    def potential_scalp_tan(self, theta):
        """Calculate sum with constants and legendres
        Parameters
        __________
        theta : ndarray [radians]
                Array of angles between electrode location and
                dipole location vectors.
        p :     ndarray [1Ee-15 Am]
                Array containing the current dipole moment for all timesteps in
                the x-, y- and z-direction.
        Returns
        _______
        pot_sum : float
                Sum factor containing brain constants, dipole and
                electrode locations and legendre polynomials.
        """
        n = np.arange(1,100)
        c4n = self.calc_c4n(n)
        d4n = self.calc_d4n(n, c4n)
        consts = c4n*(self.r/self.r4)**n + d4n*(self.r4/self.r)**(n+1)
        pot_sum = np.sum([c*lpmv(1, i, np.cos(theta)) for c,i in zip(consts,n)])
        return pot_sum

    def calc_vn(self, n):
        r_const = (self.r34 ** n - self.r43 ** (n + 1)) / ((n + 1) / n * self.r34 ** n + self.r43 ** (n + 1))
        if self.source == 'nunez' or self.source == 'srinivasan':
            v = ((n + 1) / n * self.sigma34 - r_const) / (self.sigma34 + r_const)
        else:
            v = (n / (n + 1) * self.sigma34 - r_const) / (self.sigma34 + r_const)
        return v

    def calc_yn(self, n):
        vn = self.calc_vn(n)
        r_const = (n / (n + 1) * self.r23 ** n - vn * self.r32 ** (n + 1)) / (self.r23 ** n + vn * self.r32 ** (n + 1))
        if self.source == 'nunez' or self.source == 'srinivasan':
            y = (n / (n + 1) * self.sigma23 - self.r23 ** n*(n / (n + 1)  -
                 vn * self.r32 ** (n + 1)) / (self.r23 ** n +
                 vn * self.r32 ** (n + 1))) / (self.sigma23 + r_const)
        else:
            y = (n / (n + 1) * self.sigma23 - r_const) / (self.sigma23 + r_const)
        return y

    def calc_zn(self, n):
        yn = self.calc_yn(n)
        z = (self.r12 ** n - (n + 1) / n * yn * self.r21 ** (n + 1)) / (self.r12 ** n + yn * self.r21 ** (n + 1))
        return z

    def calc_c1n(self, n):
        zn = self.calc_zn(n)
        c = ((n + 1) / n * self.sigma12 + zn) / (self.sigma12 - zn) * self.rz1**(n+1)
        return c

    def calc_c2n(self, n):
        yn = self.calc_yn(n)
        c1 = self.calc_c1n(n)
        c2 = (c1 + self.rz1**(n+1)) / (self.r12 ** n + yn * self.r21 ** (n + 1))
        return c2

    def calc_d2n(self, n, c2):
        yn = self.calc_yn(n)
        d2 = yn * c2
        return d2

    def calc_c3n(self, n):
        vn = self.calc_vn(n)
        c2 = self.calc_c2n(n)
        d2 = self.calc_d2n(n, c2)
        c3 = (c2 + d2) / (self.r23 ** n + vn * self.r32 ** (n + 1))
        return c3

    def calc_d3n(self, n, c3):
        vn = self.calc_vn(n)
        d3 = vn * c3
        return d3

    def calc_c4n(self, n):
        c3 = self.calc_c3n(n)
        d3 = self.calc_d3n(n, c3)
        c4 = (n + 1) / n * (c3 + d3) / ((n + 1) / n * self.r34 ** n + self.r43 ** (n + 1))
        if self.source == 'nunez':
            c4 = n/(n+1)*c4
        return c4

    def calc_d4n(self, n, c4):
        d4 = n / (n + 1) * c4
        return d4


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    rz1 = np.array([     0.,      0.,  78000.])
    P1 = np.array([0., 0., 1000])
    sigma=0.3
    sigmas = [sigma, sigma+1e-16, sigma+1e-14, sigma]
    r_el = np.array([0., 0., 90000.])

    dist = r_el - rz1
    cos_theta = np.dot(P1, dist)/(np.linalg.norm(dist)*np.linalg.norm(P1))
    phi_hom = 1./(4*np.pi*sigmas[0])*np.linalg.norm(P1)*cos_theta/np.sum(dist**2)*1e6

    rs = 10*np.linspace(1e4, 0.7*1e5, num=101)
    phi_homogeneous = phi_hom*np.ones(len(rs))
    phi_nunez = np.zeros(len(rs))
    phi_srinivasan = np.zeros(len(rs))
    phi_4s = np.zeros(len(rs))

    for i,r in enumerate(rs):
        radii = [79000, 80000, 85000, r]
        sphere_mod_nunez = CalcPotential4Sphere(radii, sigmas, r_el, rz1, source='nunez')
        sphere_mod_srinivasan = CalcPotential4Sphere(radii, sigmas, r_el, rz1, source='srinivasan')
        sphere_mod_4s = CalcPotential4Sphere(radii, sigmas, r_el, rz1, source=None)
        phi_nunez[i] = np.nan_to_num(sphere_mod_nunez.calc_potential(P1))
        phi_srinivasan[i] = np.nan_to_num(sphere_mod_srinivasan.calc_potential(P1))
        phi_4s[i] = np.nan_to_num(sphere_mod_4s.calc_potential(P1))

    fig = plt.figure(figsize=[5, 4])
    fig.subplots_adjust(bottom=0.14, left = 0.12)
    plt.plot(rs, phi_homogeneous, 'r', label='Infinite homogeneous space')
    plt.plot(rs, phi_nunez, 'g+', label='Nunez & Srinivasan (2006)')
    plt.plot(rs, phi_srinivasan, 'k*', label='Srinivasan (1998)')
    plt.plot(rs, phi_4s, 'b.', label='Present results - analytical')
    fig.axes[0].set_xticklabels([10, 20, 30, 40, 50, 60, 70])  # , 80, 90, 100])
    plt.xlabel(r'Scalp radius (cm)')
    plt.ylabel('Potential (nV)')
    plt.title(r'$\sigma_{skull} = \sigma_{brain} = \sigma_{csf} = \sigma_{scalp}$')
    plt.legend(frameon=False, bbox_to_anchor=(1, 0.9), fontsize=11)
    simplify_axes(fig.axes)
    # mark_subplots(fig.axes, letters='B', xpos=0., ypos = 1.1)
    plt.savefig(os.path.join('results', 'homogeneous_case.png'), dpi=300)
    # plt.show()
