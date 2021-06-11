/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_IO_RADIOMETRIES_IO_H
#define WHATPROT_IO_RADIOMETRIES_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "common/radiometry.h"
#include "common/sourced-data.h"

namespace whatprot {

void read_radiometries(const std::string& filename,
                       unsigned int* num_timesteps,
                       unsigned int* num_channels,
                       unsigned int* total_num_radiometries,
                       std::vector<Radiometry>* radiometries);

void read_radiometries_raw(const std::string& filename,
                           unsigned int* num_timesteps,
                           unsigned int* num_channels,
                           unsigned int* total_num_radiometries,
                           unsigned int* num_radiometries,
                           double** intensities);

void convert_radiometries_from_raw(unsigned int num_timesteps,
                                   unsigned int num_channels,
                                   unsigned int num_radiometries,
                                   double* intensities,
                                   std::vector<Radiometry>* radiometries);

void write_radiometries(
        const std::string& filename,
        unsigned int num_timesteps,
        unsigned int num_channels,
        const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                radiometries);

void convert_raw_from_radiometries(
        const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                radiometries,
        unsigned int radiometry_size,
        double** intensities);

void write_radiometries_raw(const std::string& filename,
                            unsigned int num_timesteps,
                            unsigned int num_channels,
                            unsigned int num_radiometries,
                            double* intensities);

void write_ys(const std::string& filename,
              const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                      radiometries);

void get_raw_ys(const std::vector<SourcedData<Radiometry, SourceCount<int>>>&
                        radiometries,
                int** ys);

void write_ys_raw(const std::string& filename,
                  unsigned int num_radiometries,
                  int* ys);

}  // namespace whatprot

#endif  // WHATPROT_IO_RADIOMETRIES_IO_H
