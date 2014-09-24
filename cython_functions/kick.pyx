
import cython
from cython.parallel cimport prange
import numpy
from cython.view cimport array as cvarray
cimport numpy
from libc.math cimport floor,fabs,sin
from openmp cimport omp_set_num_threads, omp_get_max_threads, omp_get_thread_num

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef kick(int n_rf, double[::1] beam_dE, \
         double[::1] beam_theta, \
         double[:] voltage, \
         double[:] harmonic,\
         double[:] phi_offset,\
         nthread = None): 
    
    cdef double[:] cvoltage = numpy.ascontiguousarray(voltage,dtype=numpy.float64)
    cdef double[:] charmonic = numpy.ascontiguousarray(harmonic,dtype=numpy.float64)
    cdef double[:] cphi_offset = numpy.ascontiguousarray(phi_offset,dtype=numpy.float64)
    cdef long size = beam_dE.shape[0]
    cdef int j, i
    if nthread is not None:
        if isinstance(nthread, int) and (nthread > 0):
            omp_set_num_threads(< int > nthread)
        else:
            nthread = omp_get_max_threads()
    else:
        nthread = omp_get_max_threads()
    
    for j in xrange(n_rf):
        with nogil:
            for i in prange(size):
                beam_dE[i] = beam_dE[i] + cvoltage[j] * sin(charmonic[j] * beam_theta[i] + cphi_offset[j])
    


