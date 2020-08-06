// Author: Matthew Beauregard Smith (UT Austin)
//
// IF USING MPI, YOU MUST INCLUDE <mpi.h> BEFORE INCLUDING THIS FILE.
#ifndef FLUOROSEQ_UTIL_TIME_H
#define FLUOROSEQ_UTIL_TIME_H

#ifdef MPI_VERSION  // if mpi.h has been included.

// Technically including mpi.h is redundant, since MPI_VERSION is included in
// mpi.h, so it must already be present, however this may make the file easier
// to handle in some code editors.
#include <mpi.h>

namespace fluoroseq {

double wall_time() {
    return MPI_Wtime();
}

}  // namespace fluoroseq

#else  // if mpi.h has not been included

#include <ctime>

namespace fluoroseq {

double wall_time() {
    return (double) clock() / (double) CLOCKS_PER_SEC;
}

}  // namespace fluoroseq

#endif  // MPI_VERSION

#endif  // FLUOROSEQ_UTIL_TIME_H