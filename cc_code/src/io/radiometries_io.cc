// Author: Matthew Beauregard Smith (UT Austin)
#include "radiometries_io.h"

#include <fstream>
#include <string>

#include "common/radiometry.h"

namespace fluoroseq {

namespace {
using std::ifstream;
using std::string;
}  // namespace

void read_radiometries(const string& filename,
                       int* num_timesteps,
                       int* num_channels,
                       int* num_radiometries,
                       Radiometry*** radiometries) {
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> *num_radiometries;
    *radiometries = new Radiometry*[*num_radiometries];
    for (int i = 0; i < *num_radiometries; i++) {
        (*radiometries)[i] = new Radiometry(*num_timesteps, *num_channels);
        for (int j = 0; j < (*num_timesteps) * (*num_channels); j++) {
            f >> (*radiometries)[i]->intensities[j];
        }
    }
    f.close();
}

}  // namespace fluoroseq
