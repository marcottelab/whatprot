// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_IO_RADIOMETRIES_IO_H
#define FLUOROSEQ_IO_RADIOMETRIES_IO_H

#include <string>
#include <vector>

#include "common/radiometry.h"
#include "common/sourced_data.h"

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

void write_radiometries(
        const std::string& filename,
        int num_timesteps,
        int num_channels,
        int total_num_groups,
        int group_size,
        const std::vector<
                SourcedData<Radiometry, SourceCount<int>>>& radiometries);

void convert_raw_from_radiometries(
        const std::vector<
                SourcedData<Radiometry, SourceCount<int>>>& radiometries,
        int radiometry_size,
        double** intensities);

#ifdef USE_MPI
void gather_radiometries(int total_num_blocks,
                         int block_size,
                         double** intensities);
#endif  // USE_MPI

void write_radiometries_raw(const std::string& filename,
                            int num_timesteps,
                            int num_channels,
                            int num_radiometries,
                            double* intensities);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_RADIOMETRIES_IO_H
