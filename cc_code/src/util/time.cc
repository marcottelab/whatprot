// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "time.h"

#ifdef USE_MPI
#include <mpi.h>
#else  // USE_MPI
#include <ctime>
#endif  // USE_MPI

namespace fluoroseq {

double wall_time() {
#ifdef USE_MPI
    return MPI_Wtime();
#else  // USE_MPI
    return (double) clock() / (double) CLOCKS_PER_SEC;
#endif  // USE_MPI
}

double wall_tick() {
#ifdef USE_MPI
    return MPI_Wtick();
#else  // USE_MPI
    return (double) 1.0 / (double) CLOCKS_PER_SEC;
#endif  // USE_MPI
}

unsigned int time_based_seed() {
#ifdef USE_MPI
    return (unsigned int) (wall_time() / wall_tick() + 0.1);
#else  // USE_MPI
    return (unsigned int) clock();
#endif  // USE_MPI
}

}  // namespace fluoroseq