from libcpp.vector cimport vector

# datatype definitions from the ripser library

cdef extern from "ripser.cpp":
	ctypedef struct ripserResults:
                # given a dimension d and an index k
                # the k-th pair of dimension d is
                # birth, death = births_and_death_by_dim[d][2*k], births_and_death_by_dim[d][2*k+1]
		vector[vector[float]] births_and_deaths_by_dim

                # given a dimension d and an index k
                # the k-th co-cyle of dimension d si
                # cocycle = cocycles_by_dim[d][k]
		vector[vector[vector[int]]] cocycles_by_dim

		int num_edges

        # compute given the distance matrix D
	ripserResults rips_dm(float* D, int N, int modulus, int dim_max, float threshold, bint do_cocycles)
	
        # compute given the sparse distance matrix D
	ripserResults rips_dm_sparse(int* I, int* J, float* V, int NEdges, int N, int modulus, int dim_max, float threshold, bint do_cocycles)
