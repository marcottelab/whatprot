// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_IO_RADIOMETRIES_IO_H
#define FLUOROSEQ_IO_RADIOMETRIES_IO_H

#include <string>

#include "common/radiometry.h"

namespace fluoroseq {

void read_radiometries(const std::string& filename,
                       int* num_timesteps,
                       int* num_channels,
                       int* num_radiometries,
                       Radiometry*** radiometries);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_RADIOMETRIES_IO_H