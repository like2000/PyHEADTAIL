'''
@author: Kevin Li
@date: 11.02.2014
'''


import h5py as hp
import numpy as np


from abc import ABCMeta, abstractmethod


class Monitor(object):

    @abstractmethod
    def dump(bunch):
        pass

class BunchMonitor(Monitor):

    def __init__(self, filename, n_steps, dictionary=None):

        self.h5file = hp.File(filename + '.h5', 'w')
        self.n_steps = n_steps
        self.i_steps = 0

        if dictionary:
            for key in dictionary:
                self.h5file.attrs[key] = dictionary[key]

        self.h5file.create_group('Bunch')

    # def __del__(self):

    #     self.h5file.close()
    #     print "Closed!"

    def dump(self, bunch):

        if not self.i_steps:
            n_steps = self.n_steps
            n_slices = bunch.slices.n_slices
            self.create_data(self.h5file['Bunch'], (n_steps,))
            self.write_data(bunch, np.s_[-2], self.h5file['Bunch'], self.i_steps)
        else:
            self.write_data(bunch, np.s_[-2], self.h5file['Bunch'], self.i_steps)

        self.i_steps += 1

    def create_data(self, h5group, dims):

        h5group.create_dataset("mean_x", dims)
        h5group.create_dataset("mean_xp", dims)
        h5group.create_dataset("mean_y", dims)
        h5group.create_dataset("mean_yp", dims)
        h5group.create_dataset("mean_dz", dims)
        h5group.create_dataset("mean_dp", dims)
        h5group.create_dataset("sigma_x", dims)
        h5group.create_dataset("sigma_y", dims)
        h5group.create_dataset("sigma_dz", dims)
        h5group.create_dataset("sigma_dp", dims)
        h5group.create_dataset("epsn_x", dims)
        h5group.create_dataset("epsn_y", dims)
        h5group.create_dataset("epsn_z", dims)
        h5group.create_dataset("n_macroparticles", dims)

    def write_data(self, bunch, indices, h5group, i_steps, rank=1):

        if rank == 1:
            h5group["mean_x"][i_steps] = bunch.slices.mean_x[indices]
            h5group["mean_xp"][i_steps] = bunch.slices.mean_xp[indices]
            h5group["mean_y"][i_steps] = bunch.slices.mean_y[indices]
            h5group["mean_yp"][i_steps] = bunch.slices.mean_yp[indices]
            h5group["mean_dz"][i_steps] = bunch.slices.mean_dz[indices]
            h5group["mean_dp"][i_steps] = bunch.slices.mean_dp[indices]
            h5group["sigma_x"][i_steps] = bunch.slices.sigma_x[indices]
            h5group["sigma_y"][i_steps] = bunch.slices.sigma_y[indices]
            h5group["sigma_dz"][i_steps] = bunch.slices.sigma_dz[indices]
            h5group["sigma_dp"][i_steps] = bunch.slices.sigma_dp[indices]
            h5group["epsn_x"][i_steps] = bunch.slices.epsn_x[indices]
            h5group["epsn_y"][i_steps] = bunch.slices.epsn_y[indices]
            h5group["epsn_z"][i_steps] = bunch.slices.epsn_z[indices]
            h5group["n_macroparticles"][i_steps] = bunch.slices.n_macroparticles[indices]
        elif rank == 2:
            h5group["mean_x"][:,i_steps] = bunch.slices.mean_x[indices]
            h5group["mean_xp"][:,i_steps] = bunch.slices.mean_xp[indices]
            h5group["mean_y"][:,i_steps] = bunch.slices.mean_y[indices]
            h5group["mean_yp"][:,i_steps] = bunch.slices.mean_yp[indices]
            h5group["mean_dz"][:,i_steps] = bunch.slices.mean_dz[indices]
            h5group["mean_dp"][:,i_steps] = bunch.slices.mean_dp[indices]
            h5group["sigma_x"][:,i_steps] = bunch.slices.sigma_x[indices]
            h5group["sigma_y"][:,i_steps] = bunch.slices.sigma_y[indices]
            h5group["sigma_dz"][:,i_steps] = bunch.slices.sigma_dz[indices]
            h5group["sigma_dp"][:,i_steps] = bunch.slices.sigma_dp[indices]
            h5group["epsn_x"][:,i_steps] = bunch.slices.epsn_x[indices]
            h5group["epsn_y"][:,i_steps] = bunch.slices.epsn_y[indices]
            h5group["epsn_z"][:,i_steps] = bunch.slices.epsn_z[indices]
            h5group["n_macroparticles"][:,i_steps] = bunch.slices.n_macroparticles[indices]
        else:
            raise ValueError("Rank > 2 not supported!")

class SliceMonitor(Monitor):

    def __init__(self, filename, n_steps, dictionary=None):

        self.h5file = hp.File(filename + '.h5', 'w')
        self.n_steps = n_steps
        self.i_steps = 0

        self.h5file.create_group('Bunch')
        self.h5file.create_group('LSlice')
        self.h5file.create_group('RSlice')
        self.h5file.create_group('Slices')

    # def __del__(self):

    #     self.h5file.close()
    #     print "Closed!"

    def dump(self, bunch):

        if not self.i_steps:
            n_steps = self.n_steps
            n_slices = bunch.slices.n_slices

            self.create_data(self.h5file['Bunch'], (n_steps,))
            self.create_data(self.h5file['LSlice'], (n_steps,))
            self.create_data(self.h5file['RSlice'], (n_steps,))
            self.create_data(self.h5file['Slices'], (n_slices, n_steps))

            self.write_data(bunch, np.s_[-2], self.h5file['Bunch'], self.i_steps)
            self.write_data(bunch, np.s_[0], self.h5file['LSlice'], self.i_steps)
            self.write_data(bunch, np.s_[-3], self.h5file['RSlice'], self.i_steps)
            self.write_data(bunch, np.s_[1:-3], self.h5file['Slices'], self.i_steps, rank=2)
        else:
            self.write_data(bunch, np.s_[-2], self.h5file['Bunch'], self.i_steps)
            self.write_data(bunch, np.s_[0], self.h5file['LSlice'], self.i_steps)
            self.write_data(bunch, np.s_[-3], self.h5file['RSlice'], self.i_steps)
            self.write_data(bunch, np.s_[1:-3], self.h5file['Slices'], self.i_steps, rank=2)

        self.i_steps += 1

    def create_data(self, h5group, dims):

        h5group.create_dataset("mean_x", dims)
        h5group.create_dataset("mean_xp", dims)
        h5group.create_dataset("mean_y", dims)
        h5group.create_dataset("mean_yp", dims, compression="gzip", compression_opts=9)
        h5group.create_dataset("mean_dz", dims, compression="gzip", compression_opts=9)
        h5group.create_dataset("mean_dp", dims, compression="gzip", compression_opts=9)
        h5group.create_dataset("sigma_x", dims, compression="gzip", compression_opts=9)
        h5group.create_dataset("sigma_y", dims, compression="gzip", compression_opts=9)
        h5group.create_dataset("sigma_dz", dims)
        h5group.create_dataset("sigma_dp", dims)
        h5group.create_dataset("epsn_x", dims)
        h5group.create_dataset("epsn_y", dims)
        h5group.create_dataset("epsn_z", dims)
        h5group.create_dataset("n_macroparticles", dims)

    def write_data(self, bunch, indices, h5group, i_steps, rank=1):

        if rank == 1:
            h5group["mean_x"][i_steps] = bunch.slices.mean_x[indices]
            h5group["mean_xp"][i_steps] = bunch.slices.mean_xp[indices]
            h5group["mean_y"][i_steps] = bunch.slices.mean_y[indices]
            h5group["mean_yp"][i_steps] = bunch.slices.mean_yp[indices]
            h5group["mean_dz"][i_steps] = bunch.slices.mean_dz[indices]
            h5group["mean_dp"][i_steps] = bunch.slices.mean_dp[indices]
            h5group["sigma_x"][i_steps] = bunch.slices.sigma_x[indices]
            h5group["sigma_y"][i_steps] = bunch.slices.sigma_y[indices]
            h5group["sigma_dz"][i_steps] = bunch.slices.sigma_dz[indices]
            h5group["sigma_dp"][i_steps] = bunch.slices.sigma_dp[indices]
            h5group["epsn_x"][i_steps] = bunch.slices.epsn_x[indices]
            h5group["epsn_y"][i_steps] = bunch.slices.epsn_y[indices]
            h5group["epsn_z"][i_steps] = bunch.slices.epsn_z[indices]
            h5group["n_macroparticles"][i_steps] = bunch.slices.n_macroparticles[indices]
        elif rank == 2:
            h5group["mean_x"][:,i_steps] = bunch.slices.mean_x[indices]
            h5group["mean_xp"][:,i_steps] = bunch.slices.mean_xp[indices]
            h5group["mean_y"][:,i_steps] = bunch.slices.mean_y[indices]
            h5group["mean_yp"][:,i_steps] = bunch.slices.mean_yp[indices]
            h5group["mean_dz"][:,i_steps] = bunch.slices.mean_dz[indices]
            h5group["mean_dp"][:,i_steps] = bunch.slices.mean_dp[indices]
            h5group["sigma_x"][:,i_steps] = bunch.slices.sigma_x[indices]
            h5group["sigma_y"][:,i_steps] = bunch.slices.sigma_y[indices]
            h5group["sigma_dz"][:,i_steps] = bunch.slices.sigma_dz[indices]
            h5group["sigma_dp"][:,i_steps] = bunch.slices.sigma_dp[indices]
            h5group["epsn_x"][:,i_steps] = bunch.slices.epsn_x[indices]
            h5group["epsn_y"][:,i_steps] = bunch.slices.epsn_y[indices]
            h5group["epsn_z"][:,i_steps] = bunch.slices.epsn_z[indices]
            h5group["n_macroparticles"][:,i_steps] = bunch.slices.n_macroparticles[indices]
        else:
            raise ValueError("Rank > 2 not supported!")

class ParticleMonitor(Monitor):

    def __init__(self, filename, stride=1, dictionary=None, n_steps=None):

        self.h5file = hp.File(filename + '.h5part', 'w')
        self.n_steps = n_steps
        self.i_steps = 0
        self.stride = stride

        if dictionary:
            for key in dictionary:
                self.h5file.attrs[key] = dictionary[key]
        
    def dump(self, bunch):

        if self.n_steps:
            if not self.i_steps:
                resorting_indices = np.argsort(bunch.id)[::self.stride]
                self.z0 = np.copy(bunch.dz[resorting_indices])
                n_macroparticles = len(self.z0)

                self.create_data(self.h5file, (n_macroparticles, self.n_steps))

            self.write_data(bunch, self.i_steps)
        else:
            if not self.i_steps:
                resorting_indices = np.argsort(bunch.id)[::self.stride]
                self.z0 = np.copy(bunch.dz[resorting_indices])
            n_macroparticles = len(self.z0)

            h5group = self.h5file.create_group("Step#" + str(self.i_steps))
            self.create_data(h5group, (n_macroparticles,))
            self.write_data(bunch, h5group)

        self.i_steps += 1

    def create_data(self, h5group, dims):

        if self.n_steps:
            h5group.create_dataset("x", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("xp", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("y", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("yp", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("dz", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("dp", dims, compression="gzip", compression_opts=9)
        else:
            h5group.create_dataset("x_0", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("v_0", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("x_1", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("v_1", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("x_2", dims, compression="gzip", compression_opts=9)
            h5group.create_dataset("v_2", dims, compression="gzip", compression_opts=9)

        # Do we need/want this here?
        h5group.create_dataset("id", dims, dtype=np.int, compression="gzip", compression_opts=9)
        h5group.create_dataset("c", dims, compression="gzip", compression_opts=9)

    def write_data(self, bunch, h5group):

        resorting_indices = np.argsort(bunch.id)[::self.stride]

        if self.n_steps:
            self.h5file['x'][:,h5group] = bunch.x.take(resorting_indices)
            self.h5file['xp'][:,h5group] = bunch.xp.take(resorting_indices)
            self.h5file['y'][:,h5group] = bunch.y.take(resorting_indices)
            self.h5file['yp'][:,h5group] = bunch.yp.take(resorting_indices)
            self.h5file['dz'][:,h5group] = bunch.dz.take(resorting_indices)
            self.h5file['dp'][:,h5group] = bunch.dp.take(resorting_indices)

            self.h5file['id'][:,h5group] = bunch.id.take(resorting_indices)
            self.h5file['c'][:,h5group] = self.z0
        else:
            h5group["x_0"][:] = bunch.x[resorting_indices]
            h5group["v_0"][:] = bunch.yp[resorting_indices]
            h5group["x_1"][:] = bunch.y[resorting_indices]
            h5group["v_1"][:] = bunch.yp[resorting_indices]
            h5group["x_2"][:] = bunch.dz[resorting_indices]
            h5group["v_2"][:] = bunch.dp[resorting_indices]

            # Do we need/want this here?
            h5group["id"][:] = bunch.id[resorting_indices]
            h5group["c"][:] = self.z0
