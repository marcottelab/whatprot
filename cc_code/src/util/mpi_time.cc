// Author: Matthew Beauregard Smith (UT Austin)
#include "time.h"

#include <mpi.h>

namespace fluoroseq {

double wall_time() {
    return MPI_Wtime();
}

}  // namespace fluoroseq
