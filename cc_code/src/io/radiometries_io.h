/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_IO_RADIOMETRIES_IO_H
#define FLUOROSEQ_IO_RADIOMETRIES_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
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
        const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                radiometries);

void convert_raw_from_radiometries(
        const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                radiometries,
        int radiometry_size,
        double** intensities);

#ifdef USE_MPI
void gather_radiometries(int num_radiometries,
                         int radiometry_size,
                         double** intensities,
                         int* total_num_radiometries);
#endif  // USE_MPI

void write_radiometries_raw(const std::string& filename,
                            int num_timesteps,
                            int num_channels,
                            int num_radiometries,
                            double* intensities);

void write_ys(const std::string& filename,
              const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                      radiometries);

void get_raw_ys(const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                        radiometries,
                int** ys);

#ifdef USE_MPI
void gather_ys(int num_radiometries, int** ys, int* total_num_radiometries);
#endif  // USE_MPI

void write_ys_raw(const std::string& filename, int num_radiometries, int* ys);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_RADIOMETRIES_IO_H
