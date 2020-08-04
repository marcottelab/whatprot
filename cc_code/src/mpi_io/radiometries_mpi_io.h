// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_MPI_IO_RADIOMETRIES_MPI_IO_H
#define FLUOROSEQ_MPI_IO_RADIOMETRIES_MPI_IO_H

#include <string>

#include "common/radiometry.h"

namespace fluoroseq {

void mpi_read_radiometries(const std::string& filename,
                           int* num_timesteps,
                           int* num_channels,
                           int* total_num_radiometries,
                           int* num_radiometries,
                           Radiometry*** radiometries);

void mpi_read_radiometries_master(const std::string& filename,
                                  int* num_timesteps,
                                  int* num_channels,
                                  int* total_num_radiometries,
                                  int* num_radiometries,
                                  Radiometry*** radiometries);

void mpi_read_radiometries_slave(const std::string& filename,
                                 int* num_timesteps,
                                 int* num_channels,
                                 int* total_num_radiometries,
                                 int* num_radiometries,
                                 Radiometry*** radiometries);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_MPI_IO_RADIOMETRIES_MPI_IO_H
