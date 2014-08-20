#!/usr/bin/python

from trackers.longitudinal_tracker import *
from beams.beams import *
import matplotlib as mpl

n_turns = 10001
above_transition = False

if above_transition:
	phi1 = 0
	gamma = 4 + 2.49
else:
	phi1 = np.pi
	gamma = 2.49

beta = np.sqrt(1 - gamma**-2)
a = RFSystems(np.pi * 50, [1, 2], [8000, 0*6000], [phi1, np.pi],
              [4.05**-2], m_p * beta * gamma * c * 5e-6)
beam = Particles.as_gaussian(1000, e, gamma, 1.e11, m_p, 0, 54.6408, 2, 0,
                             54.5054, 2, 15, 0.001)
n = 1
x = np.arange(n * -100, n * 100, n * 0.1)
y = np.arange(n * -0.005, n * 0.005, n * 0.00001)
X, Y = np.meshgrid(x,y)
Z = a.hamiltonian(X, Y, beam)

extent=[x.min(), x.max(), y.min(), y.max()]

import sys
import matplotlib.pyplot as plt
plt.ion()
plt.figure()

for i in xrange(n_turns):
	a.track(beam)
	if i%50 is not 0:
		continue
	sys.stdout.write("\r %d" % i)
	sys.stdout.flush()
	plt.clf()
	plt.gca()
	# im = plt.imshow(Z, interpolation="bilinear", origin="lower",
	#                  vmin=Z.min(), vmax=Z.max(), aspect="auto",
	#                  extent=[x.min(), x.max(), y.min(), y.max()])
	CS = plt.contour(X, Y, Z, 25)#, colors="black")
	plt.clabel(CS, fontsize=9, inline=1)
	# ZZ, DP = np.meshgrid(beam.z, beam.dp)
	HAM = a.hamiltonian(beam.z, beam.dp, beam)
	plt.scatter(beam.z, beam.dp, alpha=0.8, linewidth=0.1, c=HAM, marker=".")
	plt.xlim(x.min(), x.max())
	plt.ylim(y.min(), y.max())
	plt.draw()

raw_input("\nPress Enter to continue...")
