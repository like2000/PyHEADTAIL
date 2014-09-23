
import cython
from cython.parallel cimport prange
import numpy
from cython.view cimport array as cvarray
cimport numpy
import sys

from libc.math cimport floor,fabs
from openmp cimport omp_set_num_threads, omp_get_max_threads, omp_get_thread_num

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def histogram(numpy.ndarray pos not None, \
              long bins=100, \
              bin_range=None, \
              nthread=None): 
    
    cdef long  size = pos.size
    cdef double[:] cpos = numpy.ascontiguousarray(pos,dtype=numpy.float64)
    cdef numpy.ndarray[numpy.float64_t, ndim = 1] out_count = numpy.zeros(bins, dtype="float64")
    cdef double bin_edge_min, bin_edge_max
    bin_edge_min = bin_range[0]
    bin_edge_max = bin_range[1] 
    cdef double bin_width = (bin_edge_max - bin_edge_min) / (< double > (bins))
    cdef double inv_bin_width = (< double > (bins)) / (bin_edge_max - bin_edge_min)
    cdef double a = 0.0
    cdef double fbin = 0.0
    cdef double ffbin = 0.0
    cdef double tmp_count = 0.0
    cdef long   bin = 0, thread
    cdef long   i
    if nthread is not None:
        if isinstance(nthread, int) and (nthread > 0):
            omp_set_num_threads(< int > nthread)
        else:
            nthread = omp_get_max_threads()
    else:
        nthread = omp_get_max_threads()
    cdef long[:,:] big_count = cvarray(shape=(nthread, bins), itemsize=sizeof(long), format="l")
    
    big_count[:,:] = 0
    
    with nogil:
        for i in prange(size):
            a = cpos[i]
            if (a < bin_edge_min) or (a > bin_edge_max):
                continue
            fbin = (a - bin_edge_min) * inv_bin_width
            ffbin = floor(fbin)
            bin = < long > ffbin
            thread = omp_get_thread_num()
            big_count[thread, bin] += 1
        for bin in prange(bins):
            tmp_count = 0
            for thread in range(omp_get_max_threads()):
                tmp_count = tmp_count + big_count[thread, bin]
            out_count[bin] += tmp_count
            
    return  out_count


