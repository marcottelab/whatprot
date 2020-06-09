// Author: Matthew Beauregard Smith (UT Austin)
//
// This file corresponds to two possible .cc files. If not using MPI, build with
// time.cc. If using MPI, build with mpi_time.cc.
#ifndef FLUOROSEQ_UTIL_TIME_H
#define FLUOROSEQ_UTIL_TIME_H

namespace fluoroseq {

double wall_time();

}  // namespace fluoroseq

#endif  // FLUOROSEQ_UTIL_TIME_H