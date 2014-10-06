'''
@file matching
@author Kevin Li, Adrian Oeftiger
@date February 2014
@brief Module for matching transverse and longitudinal distributions
@copyright CERN
'''


from __future__ import division

from abc import ABCMeta, abstractmethod

import numpy as np
from numpy.random import RandomState

from scipy.constants import c, e
from scipy.optimize import brentq, brenth, bisect, newton
from scipy.interpolate import interp2d
from scipy.integrate import quad, fixed_quad, dblquad, cumtrapz, romb


class PhaseSpace(object):
    """Knows how to distribute particle coordinates for a beam
    according to certain distribution functions.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def generate(self, beam):
        """Creates the beam macroparticles according to a
        distribution function (depends on the implementing class).
        """
        pass


class GaussianX(PhaseSpace):
    """Horizontal Gaussian particle phase space distribution."""

    def __init__(self, sigma_x, sigma_xp, generator_seed=None):
        """Initiates the horizontal beam coordinates
        to the given Gaussian shape.
        """
        self.sigma_x  = sigma_x
        self.sigma_xp = sigma_xp

        self.random_state = RandomState()
        self.random_state.seed(generator_seed)

    @classmethod
    def from_optics(cls, alpha_x, beta_x, epsn_x, betagamma, generator_seed=None):
        """Initialise GaussianX from the given optics functions.
        beta_x is given in meters and epsn_x in micrometers.
        """

        sigma_x  = np.sqrt(beta_x * epsn_x * 1e-6 / betagamma)
        sigma_xp = sigma_x / beta_x

        return cls(sigma_x, sigma_xp, generator_seed)

    def generate(self, beam):
        beam.x = self.sigma_x * self.random_state.randn(beam.n_macroparticles)
        beam.xp = self.sigma_xp * self.random_state.randn(beam.n_macroparticles)


class GaussianY(PhaseSpace):
    """Vertical Gaussian particle phase space distribution."""

    def __init__(self, sigma_y, sigma_yp, generator_seed=None):
        """Initiates the vertical beam coordinates
        to the given Gaussian shape.
        """
        self.sigma_y  = sigma_y
        self.sigma_yp = sigma_yp

        self.random_state = RandomState()
        self.random_state.seed(generator_seed)

    @classmethod
    def from_optics(cls, alpha_y, beta_y, epsn_y, betagamma, generator_seed=None):
        """Initialise GaussianY from the given optics functions.
        beta_y is given in meters and epsn_y in micrometers.
        """

        sigma_y  = np.sqrt(beta_y * epsn_y * 1e-6 / betagamma)
        sigma_yp = sigma_y / beta_y

        return cls(sigma_y, sigma_yp, generator_seed)

    def generate(self, beam):
        beam.y = self.sigma_y * self.random_state.randn(beam.n_macroparticles)
        beam.yp = self.sigma_yp * self.random_state.randn(beam.n_macroparticles)


class GaussianZ(PhaseSpace):
    """Longitudinal Gaussian particle phase space distribution."""

    def __init__(self, sigma_z, sigma_dp, is_accepted=None, generator_seed=None):
        """Initiates the longitudinal beam coordinates to a given
        Gaussian shape. If the argument is_accepted is set to
        the is_in_separatrix(z, dp, beam) method of a RFSystems
        object (or similar), macroparticles will be initialised
        until is_accepted returns True.
        """
        self.sigma_z = sigma_z
        self.sigma_dp = sigma_dp
        self.is_accepted = is_accepted

        self.random_state = RandomState()
        self.random_state.seed(generator_seed)

    @classmethod
    def from_optics(cls, beta_z, epsn_z, p0, is_accepted=None,
                    generator_seed=None):
        """Initialise GaussianZ from the given optics functions.
        For the argument is_accepted see __init__.
        """

        sigma_z = np.sqrt(beta_z*epsn_z/(4*np.pi) * e/p0)
        sigma_dp = sigma_z / beta_z

        return cls(sigma_z, sigma_dp, is_accepted, generator_seed)

    def generate(self, beam):
        beam.z = self.sigma_z * self.random_state.randn(beam.n_macroparticles)
        beam.dp = self.sigma_dp * self.random_state.randn(beam.n_macroparticles)
        if self.is_accepted:
            self._redistribute(beam)

    def _redistribute(self, beam):
        n = beam.n_macroparticles
        z = beam.z.copy()
        dp = beam.dp.copy()

        mask_out = ~self.is_accepted(z, dp)
        while mask_out.any():
            n_gen = np.sum(mask_out)
            z[mask_out] = self.sigma_z * self.random_state.randn(n_gen)
            dp[mask_out] = self.sigma_dp * self.random_state.randn(n_gen)
            mask_out = ~self.is_accepted(z, dp)
            print 'Reiterate on non-accepted particles'

        # for i in xrange(n):
        #     while not self.is_accepted(z[i], dp[i]):
        #         z[i]  = self.sigma_z * self.random_state.randn()
        #         dp[i] = self.sigma_dp * self.random_state.randn()

        beam.z = z
        beam.dp = dp


class RFBucketMatcher(PhaseSpace):

    def __init__(self, psi, rfbucket, sigma_z=None, epsn_z=None):

        self.psi = psi
        self.H = rfbucket
        self.sigma_z = sigma_z

        self.psi_object = psi(rfbucket.hamiltonian, rfbucket.Hmax)
        self.psi = self.psi_object.function
        self.p_limits = rfbucket.separatrix

        self._compute_std = self._compute_std_cumtrapz

        if sigma_z and not epsn_z:
            self.variable = sigma_z
            self.psi_for_variable = self.psi_for_bunchlength_newton_method
        elif not sigma_z and epsn_z:
            self.variable = epsn_z
            self.psi_for_variable = self.psi_for_emittance_newton_method
        else:
            raise ValueError("Can not generate mismatched matched distribution!")

    def psi_for_emittance_newton_method(self, epsn_z):
        H = self.H

        # Maximum emittance
        self._set_psi_sigma(H.circumference)
        # zc_left, zc_right = self._get_edges_for_cut(np.exp(-2**2/2.))
        # epsn_max = self._compute_zero_quad(lambda y, x: 1, H.equihamiltonian(zc_left), H.zleft, H.zright) * 2*H.p0_reference/e
        z, dp = self._generate(50000, self.psi)
        epsn_max = self._compute_emittance(z, dp) * 4*np.pi*H.p0_reference/e
        if epsn_z > epsn_max:
            print '\n*** RMS emittance larger than bucket; using full bucket emittance', epsn_max, ' [eV s].'
            epsn_z = epsn_max*0.99
        print '\n*** Maximum RMS emittance', epsn_max, 'eV s.'

        # @profile
        def get_zc_for_epsn_z(ec):
            self._set_psi_epsn(ec)
            # zc_left, zc_right = self._get_edges_for_cut(np.exp(-2**2/2.))
            # emittance = self._compute_zero_quad(lambda y, x: 1, H.equihamiltonian(zc_left), H.zleft, H.zright) * 2*H.p0_reference/e
            z, dp = self._generate(50000, self.psi)
            emittance = self._compute_emittance(z, dp) * 4*np.pi*H.p0_reference/e
            print '... distance to target emittance:', emittance-epsn_z

            return emittance-epsn_z

        try:
            ec_bar = newton(get_zc_for_epsn_z, epsn_z*1.2, tol=1e-4, maxiter=30)
        except RuntimeError:
            print '*** WARNING: failed to converge using Newton-Raphson method. Trying classic Brent method...'
            ec_bar = brentq(get_zc_for_epsn_z, epsn_z/2, 2*epsn_max)

        self._set_psi_epsn(ec_bar)
        # zc_left, zc_right = self._get_edges_for_cut(np.exp(-2**2/2.))
        # emittance = self._compute_zero_quad(lambda y, x: 1, H.equihamiltonian(zc_left), H.zleft, H.zright) * 2*H.p0_reference/e
        z, dp = self._generate(50000, self.psi)
        emittance = self._compute_emittance(z, dp) * 4*np.pi*H.p0_reference/e
        sigma = self._compute_std(self.psi, H.separatrix, H.zleft, H.zright)

        print '\n--> Emittance:', emittance
        print '--> Bunch length:', sigma
        # H.zleft_for_eps, H.zright_for_eps = zc_left, zc_right
        # H.emittance, H.sigma = emittance, sigma

    # @profile
    def psi_for_bunchlength_newton_method(self, sigma):
        H = self.H

        # Maximum bunch length
        self._set_psi_sigma(H.circumference)
        sigma_max = self._compute_std(self.psi, H.separatrix, H.zleft, H.zright)
        if sigma > sigma_max:
            print "\n*** RMS bunch larger than bucket; using full bucket rms length", sigma_max, " m."
            sigma = sigma_max*0.99
        print '\n*** Maximum RMS bunch length', sigma_max, 'm.'

        # Width for bunch length
        def get_zc_for_sigma(zc):
            self._set_psi_sigma(zc)
            length = self._compute_std(self.psi, H.separatrix, H.zleft, H.zright)
            if np.isnan(length):
                raise ValueError
            print '... distance to target bunchlength:', length-sigma

            return length-sigma

        zc_bar = newton(get_zc_for_sigma, sigma)

        self._set_psi_sigma(zc_bar)
        # zc_left, zc_right = self._get_edges_for_cut(np.exp(-2**2/2.))
        # emittance = self._compute_zero_quad(lambda y, x: 1, H.equihamiltonian(zc_left), H.zleft, H.zright) * 2*H.p0_reference/e
        z, dp = self._generate(50000, self.psi)
        emittance = self._compute_emittance(z, dp) * 4*np.pi*H.p0_reference/e
        sigma = self._compute_std(self.psi, H.separatrix, H.zleft, H.zright)

        print '--> Emittance:', emittance
        print '\n--> Bunch length:', sigma
        # H.zleft_for_eps, H.zright_for_eps = zc_left, zc_right
        # H.emittance, H.sigma = emittance, sigma

    def generate(self, macroparticlenumber, particles=None):
        '''
        Generate a 2d phase space of n_particles particles randomly distributed
        according to the particle distribution function psi within the region
        [xmin, xmax, ymin, ymax].
        '''
        # psi = self.psi_for_variable(self.variable)
        self.psi_for_variable(self.variable)

        # Bin
        i, j = 0, 0
        nx, ny = 128, 128
        xmin, xmax = self.H.zleft, self.H.zright
        ymin, ymax = -self.H.p_max(self.H.zright), self.H.p_max(self.H.zright)
        lx = (xmax - xmin)
        ly = (ymax - ymin)

        # xx = np.linspace(xmin, xmax, nx + 1)
        # yy = np.linspace(ymin, ymax, ny + 1)
        # XX, YY = np.meshgrid(xx, yy)
        # HH = self.psi(XX, YY)
        # psi_interp = interp2d(xx, yy, HH)

        # ================================================================
        # mask_out = ~self.is_accepted(z, dp)
        # while mask_out.any():
        #     n_gen = np.sum(mask_out)
        #     z[mask_out] = self.sigma_z * self.random_state.randn(n_gen)
        #     dp[mask_out] = self.sigma_dp * self.random_state.randn(n_gen)
        #     mask_out = ~self.is_accepted(z, dp)
        #     print 'Reiterate on non-accepted particles'

        # for i in xrange(n):
        #     while not self.is_accepted(z[i], dp[i]):
        #         z[i]  = self.sigma_z * self.random_state.randn()
        #         dp[i] = self.sigma_dp * self.random_state.randn()
        # ================================================================

        n_gen = macroparticlenumber
        u = xmin + lx * np.random.random(n_gen)
        v = ymin + ly * np.random.random(n_gen)
        s = np.random.random(n_gen)
        mask_out = ~(s<self.psi(u, v))
        while mask_out.any():
            n_gen = np.sum(mask_out)
            u[mask_out] = xmin + lx * np.random.random(n_gen)
            v[mask_out] = ymin + ly * np.random.random(n_gen)
            s[mask_out] = np.random.random(n_gen)
            mask_out = ~(s<self.psi(u, v))
            # print 'Reiterate on non-accepted particles.'
            # print n_gen, '\n'

        # while j < particles.n_macroparticles:
        #     u = xmin + lx * np.random.random()
        #     v = ymin + ly * np.random.random()

        #     s = np.random.random()

        #     i += 1
        #     if s < psi_interp(u, v):
        #         x[j] = u
        #         y[j] = v
        #         # TODO: check if this does not cause problems! Setter for item does not work - not implemented!
        #         # particles.dp[j] = v
        #         j += 1

        if particles:
            particles.z = u
            particles.dp = v
            # Stick auxiliary information to particles
            particles.psi = self.psi
            particles.linedensity = self.linedensity

        return u, v, self.psi, self.linedensity

    def linedensity(self, xx):
        quad_type = fixed_quad

        L = []
        try:
            L = np.array([quad_type(lambda y: self.psi(x, y), 0, self.p_limits(x))[0] for x in xx])
            # for x in xx:
            #     y = np.linspace(0, self.p_limits(x), 129)
            #     z = self.psi(x, y)
            #     L.append(cumtrapz(z, y)[-1])
        except TypeError:
            L = quad_type(lambda y: self.psi(xx, y), 0, self.p_limits(xx))[0]
            # y = np.linspace(0, self.p_limits(xx), 129)
            # z = self.psi(xx, y)
            # L.append(cumtrapz(z, y)[-1])
        L = np.array(L)

        # # L = []
        # # for x in xx:
        # #     L.append(quad(lambda y: psi(x, y), 0, p_limits(x))[0])
        # #     if abs(x) < 0.2:
        # #         plt.ion()
        # #         zz = plt.linspace(-p_limits(x), p_limits(x))
        # #         ax1.plot(zz, psi(x, zz))
        # #         ax2.plot(x, L[-1], 'o')
        # #         plt.draw()
        # #     # L.append(quad(psi, 0, p_limits(a), args=(x))[0])
        # #     # print 2*L[-1]
        # #     # L = np.array([quad(lambda y: psi(x, y), 0, p_limits(x), args=(x))[0] for x in xx])
        # # # np.savetxt('linedensity.dat', np.array([xx, 2*np.array(L)]))

        return 2*L

    def _generate(self, macroparticlenumber, psi):

        # Bin
        i, j = 0, 0
        nx, ny = 128, 128
        xmin, xmax = self.H.zleft, self.H.zright
        ymin, ymax = -self.H.p_max(self.H.zright), self.H.p_max(self.H.zright)
        lx = (xmax - xmin)
        ly = (ymax - ymin)

        n_gen = macroparticlenumber
        u = xmin + lx * np.random.random(n_gen)
        v = ymin + ly * np.random.random(n_gen)
        s = np.random.random(n_gen)
        mask_out = ~(s<self.psi(u, v))
        while mask_out.any():
            n_gen = np.sum(mask_out)
            u[mask_out] = xmin + lx * np.random.random(n_gen)
            v[mask_out] = ymin + ly * np.random.random(n_gen)
            s[mask_out] = np.random.random(n_gen)
            mask_out = ~(s<self.psi(u, v))

        return u, v

    def _set_psi_sigma(self, sigma):
        self.psi_object.H0 = self.H.H0_from_sigma(sigma)

    def _set_psi_epsn(self, epsn):
        self.psi_object.H0 = self.H.H0_from_epsn(epsn)

    # @profile
    def _get_edges_for_cut(self, h_cut):
        zz = np.linspace(self.H.zmin, self.H.zmax, 128)
        ll = self.linedensity(zz)
        lmax = np.amax(ll)
        # plt.plot(zz, linedensity(zz)/lmax, 'r', lw=2)
        # plt.plot(zz, psi(zz, 0))
        # plt.axhline(h_cut)
        # plt.axvline(zcut_bar)
        # plt.show()
        return self.H._get_zero_crossings(lambda x: self.linedensity(x) - h_cut*lmax)

    def _compute_emittance(self, z, dp):
        var_z    = np.var(z)
        var_dp   = np.var(dp)
        mean_zdp = np.mean( (z-np.mean(z)) * (dp-np.mean(dp)) )

        return np.sqrt(var_z*var_dp - mean_zdp**2)

    def _compute_zero_quad(self, psi, p_sep, xmin, xmax):
        '''
        Compute the variance of the distribution function psi from xmin to xmax
        along the contours p_sep using numerical integration methods.
        '''

        Q, error = dblquad(lambda y, x: psi(x, y), xmin, xmax,
                    lambda x: 0, lambda x: p_sep(x))

        return Q

    def _compute_mean_quad(self, psi, p_sep, xmin, xmax):
        '''
        Compute the variance of the distribution function psi from xmin to xmax
        along the contours p_sep using numerical integration methods.
        '''

        Q = self._compute_zero_quad(psi, p_sep, xmin, xmax)
        M, error = dblquad(lambda y, x: x * psi(x, y), xmin, xmax,
                    lambda x: 0, lambda x: p_sep(x))

        return M/Q

    def _compute_std_quad(self, psi, p_sep, xmin, xmax):
        '''
        Compute the variance of the distribution function psi from xmin to xmax
        along the contours p_sep using numerical integration methods.
        '''

        Q = self._compute_zero_quad(psi, p_sep, xmin, xmax)
        M = self._compute_mean_quad(psi, p_sep, xmin, xmax)
        V, error = dblquad(lambda y, x: (x-M) ** 2 * psi(x, y), xmin, xmax,
                    lambda x: 0, lambda x: p_sep(x))

        return np.sqrt(V/Q)

    def _compute_zero_cumtrapz(self, psi, p_sep, xmin, xmax):

        x_arr = np.linspace(xmin, xmax, 257)
        dx = x_arr[1] - x_arr[0]

        Q = 0
        for x in x_arr:
            y = np.linspace(0, p_sep(x), 257)
            z = psi(x, y)
            Q += cumtrapz(z, y)[-1]
        Q *= dx

        return Q

    def _compute_mean_cumtrapz(self, psi, p_sep, xmin, xmax):

        Q = self._compute_zero_cumtrapz(psi, p_sep, xmin, xmax)

        x_arr = np.linspace(xmin, xmax, 257)
        dx = x_arr[1] - x_arr[0]

        M = 0
        for x in x_arr:
            y = np.linspace(0, p_sep(x), 257)
            z = x * psi(x, y)
            M += cumtrapz(z, y)[-1]
        M *= dx

        return M/Q

    def _compute_std_cumtrapz(self, psi, p_sep, xmin, xmax):
        '''
        Compute the variance of the distribution function psi from xmin to xmax
        along the contours p_sep using numerical integration methods.
        '''

        Q = self._compute_zero_cumtrapz(psi, p_sep, xmin, xmax)
        M = self._compute_mean_cumtrapz(psi, p_sep, xmin, xmax)

        x_arr = np.linspace(xmin, xmax, 257)
        dx = x_arr[1] - x_arr[0]

        V = 0
        for x in x_arr:
            y = np.linspace(0, p_sep(x), 257)
            z = (x-M)**2 * psi(x, y)
            V += cumtrapz(z, y)[-1]
        V *= dx

        return np.sqrt(V/Q)

    def _compute_std_romberg(self, psi, p_sep, xmin, xmax):
        '''
        Compute the variance of the distribution function psi from xmin to xmax
        along the contours p_sep using numerical integration methods.
        '''

        x_arr = np.linspace(xmin, xmax, 257)
        dx = x_arr[1] - x_arr[0]

        Q, V = 0, 0
        for x in x_arr:
            y = np.linspace(0, p_sep(x), 257)
            dy = y[1] - y[0]
            z = psi(x, y)
            Q += romb(z, dy)
            z = x**2 * psi(x, y)
            V += romb(z, dy)
        Q *= dx
        V *= dx

        return np.sqrt(V/Q)


class StationaryExponential(object):

    def __init__(self, H, Hmax=None, width=1000, Hcut=0):
        self.H = H
        self.H0 = 1
        if not Hmax:
            self.Hmax = H(0, 0)
        else:
            self.Hmax = Hmax
        self.Hcut = Hcut
        self.width = width

    def function(self, z, dp):
        # psi = np.exp((self.H(z, dp)) / (self.width*self.Hmax)) - 1
        # psi_offset = np.exp(self.Hcut / (self.width*self.Hmax)) - 1
        # psi_norm = (np.exp(1/self.width) - 1) - psi_offset
        # return ( (psi-psi_offset) / psi_norm ).clip(min=0)

        # psi = np.exp( (self.H(z, dp)-self.Hcut).clip(min=0) / (self.width*self.Hmax)) - 1
        # psi_norm = np.exp( (self.Hmax-0*self.Hcut) / (self.width*self.Hmax) ) - 1
        # psi = np.exp( -self.H(z, dp).clip(min=0)/(self.width*self.Hmax) ) - 1
        # psi_norm = np.exp( -self.Hmax/(self.width*self.Hmax) ) - 1

        psi = np.exp(self.H(z, dp).clip(min=0)/self.H0) - 1
        psi_norm = np.exp(self.Hmax/self.H0) - 1
        return psi/psi_norm
