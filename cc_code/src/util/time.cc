// Author: Matthew Beauregard Smith (UT Austin)
#include "time.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <ctime>
#endif

namespace fluoroseq {

double wall_time() {
    #ifdef _OPENMP
    return omp_get_wtime();
    #else
    return (double) clock() / (double) CLOCKS_PER_SEC;
    #endif
}

}  // namespace fluoroseq
