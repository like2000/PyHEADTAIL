from __future__ import division
'''
Created on 08.01.2014

@author: Kevin Li
'''


import numpy as np


import copy
from solvers.grid import *
from solvers.compute_potential_fgreenm2m import compute_potential_fgreenm2m, compute_potential_fgreenp2m


class PoissonFFT(UniformGrid):
    '''
    FFT Poisson solver operates on a grid
    '''

    # @profile
    def __init__(self, *args, **kwargs):
        '''
        Constructor
        '''
        super(PoissonFFT, self).__init__(*args, **kwargs)

        self.tmprho = np.zeros((2 * self.ny, 2 * self.nx))
        self.fgreen = np.zeros((2 * self.ny, 2 * self.nx))

        mx = -self.dx / 2 + np.arange(self.nx + 1) * self.dx
        my = -self.dy / 2 + np.arange(self.ny + 1) * self.dy
        x, y = np.meshgrid(mx, my)
        r2 = x ** 2 + y ** 2
        # Antiderivative
        tmpfgreen = -1 / 2 * (-3 * x * y + x * y * np.log(r2)
                   + x * x * np.arctan(y / x) + y * y * np.arctan(x / y)) # * 2 / dx / dy

        # Integration and circular Green's function
        self.fgreen[:self.ny, :self.nx] = tmpfgreen[1:, 1:] + tmpfgreen[:-1, :-1] - tmpfgreen[1:, :-1] - tmpfgreen[:-1, 1:]
        self.fgreen[self.ny:, :self.nx] = self.fgreen[self.ny:0:-1, :self.nx]
        self.fgreen[:self.ny, self.nx:] = self.fgreen[:self.ny, self.nx:0:-1]
        self.fgreen[self.ny:, self.nx:] = self.fgreen[self.ny:0:-1, self.nx:0:-1]
        # # Would expect to be fully symmetric
        # self.fgreen[self.nx:, :self.ny] = self.fgreen[self.nx - 1::-1, :self.ny]
        # self.fgreen[:self.nx, self.ny:] = self.fgreen[:self.nx, self.ny - 1::-1]
        # self.fgreen[self.nx:, self.ny:] = self.fgreen[self.nx - 1::-1, self.ny - 1::-1]

        from types import MethodType
        PoissonFFT.compute_potential_fgreenm2m = MethodType(compute_potential_fgreenm2m, None, PoissonFFT)
        PoissonFFT.compute_potential_fgreenp2m = MethodType(compute_potential_fgreenp2m, None, PoissonFFT)

        try:
            import pyfftw

            # Arrays
            self.fftw_fgreen = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
            self.fftw_rho = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
            self.fftw_phi = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
            self.ifftw_fgreen = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
            self.ifftw_rho = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
            self.ifftw_phi = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')

            # Plans
            self.pfftw_fgreen = pyfftw.FFTW(self.fftw_fgreen, self.ifftw_fgreen, axes=(0,1))#, flags=('FFTW_EXHAUSTIVE',), threads=1)
            self.pfftw_rho = pyfftw.FFTW(self.fftw_rho, self.ifftw_rho, axes=(0,1))#, flags=('FFTW_EXHAUSTIVE',), threads=1)
            self.pfftw_phi = pyfftw.FFTW(self.ifftw_phi, self.fftw_phi, axes=(0,1), direction='FFTW_BACKWARD')#, flags=('FFTW_EXHAUSTIVE',), threads=1)

            self.compute_potential = self.compute_potential_fftw
        except ImportError:
            print '*** WARNING: pyfftw not available. Falling back to NumPy FFT.'
            self.compute_potential = self.compute_potential_numpy

    # @profile
    def compute_potential_numpy(self, rho, phi):

        self.tmprho[:self.ny, :self.nx] = rho

        fftphi = np.fft.fft2(self.tmprho) * np.fft.fft2(self.fgreen)

        tmpphi = np.fft.ifft2(fftphi)
        phi[:] = np.abs(tmpphi[:self.ny, :self.nx])

        # for (size_t j=0; j<np; j++)
        # {
        #     tmpphi[j] = std::sqrt(fftw_phi[j][0] * fftw_phi[j][0]
        #               + fftw_phi[j][1] * fftw_phi[j][1]);
        #     tmpphi[j] *= norm; // FFT specific
        # }

    # @profile
    def compute_potential_fftw(self, rho, phi):

        # // FFT solver
        # fftw_fgreen = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np);
        # fftw_phi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np);
        # fftw_rho = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np);
        # fftw_fgreen_T = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np);
        # fftw_phi_T = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np);
        # fftw_rho_T = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np);

        # p_fgreen = fftw_plan_dft_2d(2 * n_points_y, 2 * n_points_x,
        #                      fftw_fgreen, fftw_fgreen_T,
        #                      FFTW_FORWARD, FFTW_MEASURE);
        # p_rho = fftw_plan_dft_2d(2 * n_points_y, 2 * n_points_x,
        #                      fftw_rho, fftw_rho_T,
        #                      FFTW_FORWARD, FFTW_MEASURE);
        # p_phi_I = fftw_plan_dft_2d(2 * n_points_y, 2 * n_points_x,
        #                      fftw_phi_T, fftw_phi,
        #                      FFTW_BACKWARD, FFTW_MEASURE);

        # # Arrays
        # fftw_fgreen = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
        # fftw_rho = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
        # fftw_phi = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
        # ifftw_fgreen = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
        # ifftw_rho = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')
        # ifftw_phi = pyfftw.n_byte_align_empty((2 * self.ny, 2 * self.nx), 16, 'complex128')

        # # Plans
        # pfftw_fgreen = pyfftw.FFTW(fftw_fgreen, ifftw_fgreen, axes=(0,1))
        # pfftw_rho = pyfftw.FFTW(fftw_rho, ifftw_rho, axes=(0,1))
        # pfftw_phi = pyfftw.FFTW(ifftw_phi, fftw_phi, axes=(0,1), direction='FFTW_BACKWARD')

        # Fill arrays
        self.fftw_fgreen[:] = self.fgreen
        self.fftw_rho[:self.ny, :self.nx] = rho

        # Solve
        tmpfgreen = self.pfftw_fgreen()
        tmprho = self.pfftw_rho()

        self.ifftw_phi[:] = np.asarray(tmprho) * np.asarray(tmpfgreen)
        # self.ifftw_phi[:] = np.dot(tmprho, tmpfgreen)
        tmpphi = self.pfftw_phi()
        phi[:] = np.abs(tmpphi[:self.ny, :self.nx])

        # for (size_t j=0; j<np; j++)
        # {
        #     tmpphi[j] = std::sqrt(fftw_phi[j][0] * fftw_phi[j][0]
        #               + fftw_phi[j][1] * fftw_phi[j][1]);
        #     tmpphi[j] *= norm; // FFT specific
        # }

    def compute_fields(self, phi, ex, ey):

        ey[:], ex[:] = np.gradient(phi, self.dy, self.dx)
        ex[:] *= -1
        ey[:] *= -1

    # @profile
    def py_green_m2m(self):

        self.phi = np.zeros((self.ny, self.nx))

        da = self.dx * self.dy
        x, y = self.x.flatten(), self.y.flatten()
        phi, rho = self.phi.flatten(), self.rho.flatten()

        for i in xrange(self.n_points):
            xi, yi = x[i], y[i]
            for j in xrange(self.n_points):
                xj, yj = x[j], y[j]
                r2 = (xi - xj) ** 2 + (yi - yj) ** 2 + 1e-6
                phi[i] += (-rho[j] * da * 1 / 2 * np.log(r2)) # * self.dx * self.dy / (2 * np.pi))

        self.phi = np.reshape(phi, (self.ny, self.nx))
        #     for j in points:
        #         r = j - i + 1e-6
        #         phi += np.log(r)

        # self.phi = 0

# class PoissonGreen(RectilinearGrid2D):

#     def __init__(self, n_points_x, n_points_y, dim_x, dim_y):
#         # Inheritance
#         RectilinearGrid2D.__init__(self, "regular", n_points_x, n_points_y,
#                                    dim_x, dim_y)

#         # Allocation
#         self.dx = 2 * abs(self.dim_x) / (self.n_points_x - 1)
#         self.dy = 2 * abs(self.dim_y) / (self.n_points_y - 1)
#         self.fgreen = np.zeros((4 * self.n_points_x * self.n_points_y))
#         self.fgreeni = np.zeros((4 * self.n_points_x * self.n_points_y))

#         # Initialisation
#         self.initd()
#         self.baseGreen()
#         self.transformedGreen()
#         self.integratedGreen()

#     def initd(self):
#         ''' Double grid initialisation'''
#         mx = [self.dim_x + (i + 1) * self.dx for i in range(self.n_points_x)]
#         my = [self.dim_y + (j + 1) * self.dy for j in range(self.n_points_y)]

#         self.mx = np.append(self.mx, mx)
#         self.my = np.append(self.my, my)

#     def baseGreen(self):
#         '''
#         Standard Green's function: log(r)
#         '''
#         tmpfgreen = np.zeros((self.n_points_x, self.n_points_y))
#         for i in range(self.n_points_x):
#             for j in range(self.n_points_y):
#                 if i != 0 or j != 0:
#                     x = i * self.dx
#                     y = j * self.dy
#                     r2 = x ** 2 + y ** 2
#                     tmpfgreen[i][j] = -1 / 2. * np.log(r2)
#         tmpfgreen[0][0] = tmpfgreen[1][0] / 2. + tmpfgreen[0][1] / 2.

#         # Base region
#         for i in range(self.n_points_x):
#             for j in range(self.n_points_y):
#                 k = 2 * self.n_points_y * i + j
#                 self.fgreen[k] = tmpfgreen[i][j]
#         # Mirror x
#         for i in range(self.n_points_x):
#             for j in range(1, self.n_points_y):
#                 k = (2 * self.n_points_y) * (i + 1) - j
#                 self.fgreen[k] = tmpfgreen[i][j]
#         # Mirror y
#         for i in range(1, self.n_points_x):
#             for j in range(self.n_points_y):
#                 k = ((2 * self.n_points_y)
#                    * (2 * self.n_points_x - i) + j)
#                 self.fgreen[k] = tmpfgreen[i][j]
#         # Mirror xy
#         for i in range(1, self.n_points_x):
#             for j in range(1, self.n_points_y):
#                 k = ((2 * self.n_points_y)
#                    * (2 * self.n_points_x - (i - 1)) - j)
#                 self.fgreen[k] = tmpfgreen[i][j]

#     def transformedGreen(self):
#         '''
#         Green's function in Fourier space: 1 / (k ** 2)
#         '''
#         self.G = np.zeros((4 * self.n_points_x * self.n_points_y))
#         kx, ky = np.pi / abs(self.dim_x), np.pi / abs(self.dim_y)

#         tmpfgreen = np.zeros((self.n_points_x, self.n_points_y))
#         for i in range(self.n_points_x):
#             for j in range(self.n_points_y):
#                 x = i * kx
#                 y = j * ky
#                 tmpfgreen[i, j] = -1. / np.sqrt(x ** 2 + y ** 2 + 1e-5 ** 2)

#         # Base region
#         for i in range(self.n_points_x):
#             for j in range(self.n_points_y):
#                 k = (2 * self.n_points_y - 1) * i + j
#                 self.G[k] = tmpfgreen[i][j]
#         # Mirror x
#         for i in range(self.n_points_x - 1):
#             for j in range(self.n_points_y):
#                 k = (2 * self.n_points_y - 1) \
#                   * (2 * self.n_points_x - 2 - i) + j
#                 self.G[k] = tmpfgreen[i][j]
#         # Mirror y
#         for i in range(self.n_points_x):
#             for j in range(self.n_points_y - 1):
#                 k = (2 * self.n_points_y - 1) * (i + 1) - j - 1
#                 self.G[k] = tmpfgreen[i][j]
#         # Mirror xy
#         for i in range(self.n_points_x - 1):
#             for j in range(self.n_points_y - 1):
#                 k = (2 * self.n_points_y - 1) \
#                   * (2 * self.n_points_x - 2 - i + 1) - j - 1
#                 self.G[k] = tmpfgreen[i][j]

#     def integratedGreen(self):

#         tmpfgreen = np.zeros((self.n_points_x + 1, self.n_points_y + 1))
#         for i in range(self.n_points_x + 1):
#             for j in range(self.n_points_y + 1):
#                 x = -self.dx / 2. + i * self.dx
#                 y = -self.dy / 2. + j * self.dy
#                 r2 = x ** 2 + y ** 2
#                 tmpfgreen[i][j] = -1 / 2. * (-3 * x * y + x * y * log(r2)
#                     + x * x * np.arctan(y / x) + y * y * np.arctan(x / y))

#         # Base region
#         for i in range(self.n_points_x):
#             for j in range(self.n_points_y):
#                 k = i * 2 * self.n_points_y + j
#                 self.fgreeni[k] += tmpfgreen[i][j]
#                 self.fgreeni[k] -= tmpfgreen[i + 1][j]
#                 self.fgreeni[k] -= tmpfgreen[i][j + 1]
#                 self.fgreeni[k] += tmpfgreen[i + 1][j + 1]
#         # Mirror x
#         for i in range(self.n_points_x):
#             for j in range(1, self.n_points_y):
#                 k = (2 * self.n_points_y) * (i + 1)  - j
#                 self.fgreeni[k] += tmpfgreen[i][j]
#                 self.fgreeni[k] -= tmpfgreen[i + 1][j]
#                 self.fgreeni[k] -= tmpfgreen[i][j + 1]
#                 self.fgreeni[k] += tmpfgreen[i + 1][j + 1]
#         # Mirror y
#         for i in range(1, self.n_points_x):
#             for j in range(self.n_points_y):
#                 k = ((2 * self.n_points_x)
#                    * (2 * self.n_points_y - i) + j)
#                 self.fgreeni[k] += tmpfgreen[i][j]
#                 self.fgreeni[k] -= tmpfgreen[i + 1][j]
#                 self.fgreeni[k] -= tmpfgreen[i][j + 1]
#                 self.fgreeni[k] += tmpfgreen[i + 1][j + 1]
#         # Mirror xy
#         for i in range(1, self.n_points_x):
#             for j in range(1, self.n_points_y):
#                 k = ((2 * self.n_points_x)
#                    * (2 * self.n_points_y - (i - 1)) - (j))
#                 self.fgreeni[k] += tmpfgreen[i][j]
#                 self.fgreeni[k] -= tmpfgreen[i + 1][j]
#                 self.fgreeni[k] -= tmpfgreen[i][j + 1]
#                 self.fgreeni[k] += tmpfgreen[i + 1][j + 1]
# #        # Mirror x
# #        for i in range(self.n_points_x):
# #            for j in range(self.n_points_y - 1):
# #                k = (2 * self.n_points_y) * (i + 1)  - (j + 1)
# #                self.fgreeni[k] += tmpfgreen[i][j]
# #                self.fgreeni[k] -= tmpfgreen[i + 1][j]
# #                self.fgreeni[k] -= tmpfgreen[i][j + 1]
# #                self.fgreeni[k] += tmpfgreen[i + 1][j + 1]
# #        # Mirror y
# #        for i in range(self.n_points_x - 1):
# #            for j in range(self.n_points_y):
# #                k = ((2 * self.n_points_y)
# #                   * (2 * self.n_points_x - (i + 1)) + j)
# #                self.fgreeni[k] += tmpfgreen[i][j]
# #                self.fgreeni[k] -= tmpfgreen[i + 1][j]
# #                self.fgreeni[k] -= tmpfgreen[i][j + 1]
# #                self.fgreeni[k] += tmpfgreen[i + 1][j + 1]
# #        # Mirror xy
# #        for i in range(self.n_points_x - 1):
# #            for j in range(self.n_points_y - 1):
# #                k = ((2 * self.n_points_y)
# #                   * (2 * self.n_points_x - i) - (j + 1))
# #                self.fgreeni[k] += tmpfgreen[i][j]
# #                self.fgreeni[k] -= tmpfgreen[i + 1][j]
# #                self.fgreeni[k] -= tmpfgreen[i][j + 1]
# #                self.fgreeni[k] += tmpfgreen[i + 1][j + 1]

#     def m2mGreen(self, b):
#         '''
#         Points-to-mesh solver
#         '''
#         phi = np.zeros((self.n_points_x, self.n_points_y))
#         rho = np.reshape(b, (self.n_points_x, self.n_points_y))

#         for i in range(self.n_points_x):
#             for j in range(self.n_points_y):
#                 x0 = -abs(self.dim_x) + i * self.dx
#                 y0 = -abs(self.dim_y) + j * self.dy
#                 for k in range(self.n_points_x):
#                     for l in range(self.n_points_y):
#                         x = -abs(self.dim_x) + k * self.dx
#                         y = -abs(self.dim_y) + l * self.dy
#                         r2 = (x0 - x) ** 2 + (y0 - y) ** 2 + 1e-3 ** 2
#                         phi[i, j] += (-rho[k, l] * 1 / 2. * np.log(r2)
#                             * self.dx * self.dy / (2 * np.pi))

#         # Plot
#         figure(10)
#         XX, YY = np.meshgrid(self.mx[:self.n_points_x],
#                              self.my[:self.n_points_y])
#         XX, YY = XX.T, YY.T
#         ct = contour(XX, YY, phi, 40)
# #        clabel(ct, manual=True)
#         gca().set_aspect('equal')

#         return phi

#     def p2mGreen(self, beam):
#         '''
#         Piont-to-point solver
#         '''
#         phi = np.zeros((self.n_points_x, self.n_points_y))

#         for i in range(self.n_points_y):
#             for j in range(self.n_points_y):
#                 x0 = -abs(self.dim_x) + i * self.dx
#                 y0 = -abs(self.dim_y) + j * self.dy
#                 for k in range(int(beam.n_particles)):
#                         xp = beam.x[k]
#                         yp = beam.y[k]
#                         r2 = (x0 - xp) ** 2 + (y0 - yp) ** 2
#                         phi[i, j] += -1 / 2. * np.log(r2)

#         # Plot
#         figure(11)
#         XX, YY = np.meshgrid(self.mx[:self.n_points_x],
#                              self.my[:self.n_points_y])
#         XX, YY = XX.T, YY.T
#         contour(XX, YY, phi, 40)
# #        gca().set_xlim(np.amin(self.mx), abs(np.amin(self.mx)))
# #        gca().set_ylim(np.amin(self.my), abs(np.amin(self.my)))
#         gca().set_aspect('equal')

#         return phi

#     def solveFFT(self, b):
#         b *= 0
#         b[528] = 96100
        
#         nx = 2 * self.n_points_x
#         ny = 2 * self.n_points_y

#         rho = np.zeros((nx, ny))
#         rho[:self.n_points_x, :self.n_points_y] = np.reshape(
#             b, (self.n_points_x, self.n_points_y))

#         fgreen = np.reshape(self.fgreen, (nx, ny))
#         fgreeni = np.reshape(self.fgreeni, (nx, ny)) * 1. / (self.dx * self.dy)
#         G = np.reshape(self.G, (nx, ny))

# #         Mysterious delete
# #        fgreeni = np.insert(fgreeni, self.n_points_x, 0, 0)
# #        fgreeni = np.insert(fgreeni, self.n_points_y, 0, 1)
# #        fgreeni = np.delete(fgreeni, fgreeni.shape[0] - 1, 0)
# #        fgreeni = np.delete(fgreeni, fgreeni.shape[1] - 1, 1)

# #        fgreen = np.delete(fgreen, self.n_points_x, 0)
# #        fgreen = np.delete(fgreen, fgreen.shape[0] - 1, 0)
# #        fgreen = np.delete(fgreen, self.n_points_y, 1)
# #        fgreen = np.delete(fgreen, fgreen.shape[1] - 1, 1)
# #
# #        fgreeni = np.delete(fgreeni, self.n_points_x, 0)
# #        fgreeni = np.delete(fgreeni, fgreeni.shape[0] - 1, 0)
# #        fgreeni = np.delete(fgreeni, self.n_points_y, 1)
# #        fgreeni = np.delete(fgreeni, fgreeni.shape[1] - 1, 1)
# #
# #        G = np.delete(G, self.n_points_x, 0)
# #        G = np.delete(G, G.shape[0] - 1, 0)
# #        G = np.delete(G, self.n_points_y, 1)
# #        G = np.delete(G, G.shape[1] - 1, 1)

#         # Perform FFT
#         rhofft = np.fft.fft2(rho)
#         fgreenfft = np.fft.fft2(fgreen)
#         fgreenffti = np.fft.fft2(fgreeni)

#         phi = np.fft.ifft2(rhofft * fgreenfft)#[:self.n_points_x, :self.n_points_y]
#         phii = np.fft.ifft2(rhofft * fgreenffti)#[:self.n_points_x, :self.n_points_y]
#         phiG = np.fft.ifft2(rhofft * G)#[:self.n_points_x, :self.n_points_y]
        
#         # Rescale for volume correction
#         phi *= self.dx * self.dy / (2 * np.pi)
#         phii *= self.dx * self.dy / (2 * np.pi)
#         phiG *= self.dx * self.dy / (2 * np.pi)

# #        Green = np.fft.ifft2(G)

#         # Plot
#         mx = self.mx
#         my = self.my
# #        mx = self.mx[:self.n_points_x]
# #        my = self.my[:self.n_points_y]
# #        [axvline(i, c='r', lw=0.5) for i in mx]
# #        [axhline(i, c='r', lw=0.5) for i in my]
#         XX, YY = np.meshgrid(mx, my)
#         XX, YY = XX.T, YY.T

# #        figure(1)
# #        ct = contour(XX, YY, abs(phi), 40)
# #        scatter(XX, YY, c=rho, s=0.001 * rho, lw=0)
# #        gca().set_xlim(np.amin(mx), abs(np.amin(mx)))
# #        gca().set_ylim(np.amin(my), abs(np.amin(my)))
# #        gca().set_aspect('equal')
# #        figure(2)
#         contour(XX, YY, abs(phii), 40)
# #        scatter(XX, YY, c=rho, s=0.001 * rho, lw=0)
# #        gca().set_xlim(np.amin(mx), abs(np.amin(mx)))
# #        gca().set_ylim(np.amin(my), abs(np.amin(my)))
# #        gca().set_aspect('equal')
# #        figure(3)
# #        contour(XX, YY, abs(fgreeni), 40)
# #        gca().set_xlim(np.amin(mx), abs(np.amin(mx)))
# #        gca().set_ylim(np.amin(my), abs(np.amin(my)))
# #        gca().set_aspect('equal')

# #        savetxt("myfgreen.dat", fgreeni)

#         return phii
