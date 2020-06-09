// Author: Matthew Beauregard Smith (UT Austin)
#include "time.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <mpi.h>
#endif

namespace fluoroseq {

double wall_time() {
    #ifdef _OPENMP
    return omp_get_wtime();
    #else
    return MPI_Wtime();
    #endif
}

}  // namespace fluoroseq
