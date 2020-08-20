// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_IO_RADIOMETRIES_IO_H
#define FLUOROSEQ_IO_RADIOMETRIES_IO_H

#include <string>
#include <vector>

#include "common/radiometry.h"

namespace fluoroseq {

void read_radiometries(const std::string& filename,
                       int* num_timesteps,
                       int* num_channels,
                       int* total_num_radiometries,
                       std::vector<Radiometry>* radiometries);

void read_radiometries_raw(const std::string& filename,
                           int* num_timesteps,
                           int* num_channels,
                           int* total_num_radiometries,
                           int* num_radiometries,
                           double** intensities);

#ifdef USE_MPI
void scatter_radiometries(int* num_timesteps,
                          int* num_channels,
                          int* total_num_radiometries,
                          int* num_radiometries,
                          double** intensities);
#endif  // USE_MPI

void convert_radiometries_from_raw(int num_timesteps,
                                   int num_channels,
                                   int num_radiometries,
                                   double* intensities,
                                   std::vector<Radiometry>* radiometries);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_RADIOMETRIES_IO_H
