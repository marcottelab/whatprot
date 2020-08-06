// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI before including this header.
#ifndef FLUOROSEQ_UTIL_TIME_H
#define FLUOROSEQ_UTIL_TIME_H

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

}  // namespace fluoroseq

#endif  // FLUOROSEQ_UTIL_TIME_H