'''
**Module to set up MPI parallelization used in the longitudinal tracker.**

:Authors: **Helga Timko**
'''

from mpi4py import MPI


class MPI_Config(object):
    '''
    *Configuring MPI*
    '''
        
    def __init__(self):
        
        try:   
            self.mpi_comm = MPI.COMM_WORLD
            self.mpi_size = self.mpi_comm.Get_rank()
            self.mpi_rank = self.mpi_comm.Get_rank()
        except:
            raise RuntimeError('ERROR: Number of cores for longitudinal tracker not recognized!')

