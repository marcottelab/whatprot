// Author: Matthew Beauregard Smith (UT Austin)
//
// IF USING MPI, YOU MUST INCLUDE <mpi.h> BEFORE INCLUDING THIS FILE.
#ifndef FLUOROSEQ_UTIL_TIME_H
#define FLUOROSEQ_UTIL_TIME_H

// Technically including mpi.h again for the "mpi.h has been included" case is
// redundant, since MPI_VERSION is included in mpi.h, so it must already be
// present, however this may make the file easier to handle in some code
// editors.
#ifdef MPI_VERSION  // if mpi.h has been included.
#include <mpi.h>
#else  // if mpi.h has not been included
#include <ctime>
#endif  // MPI_VERSION

namespace fluoroseq {

double wall_time() {
#ifdef MPI_VERSION  // if mpi.h has been included.
    return MPI_Wtime();
#else  // if mpi.h has not been included.
    return (double) clock() / (double) CLOCKS_PER_SEC;
#endif  // MPI_VERSION
}

}  // namespace fluoroseq

#endif  // FLUOROSEQ_UTIL_TIME_H