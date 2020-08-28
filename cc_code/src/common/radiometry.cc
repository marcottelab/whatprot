// Author: Matthew Beauregard Smith
#include "radiometry.h"

#include <algorithm>

namespace fluoroseq {

namespace {
using std::copy;
}

Radiometry::Radiometry(int num_timesteps, int num_channels)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    intensities = new double[num_timesteps * num_channels];
}

Radiometry::Radiometry(const Radiometry& other)
        : num_timesteps(other.num_timesteps),
          num_channels(other.num_channels) {
    intensities = new double[num_timesteps * num_channels];
    copy(other.intensities,
         &other.intensities[num_timesteps * num_channels],
         intensities);
}

Radiometry::~Radiometry() {
    delete[] intensities;
}

double& Radiometry::operator()(int t, int c) {
    return intensities[t * num_channels + c];
}

double Radiometry::operator()(int t, int c) const {
    return intensities[t * num_channels + c];
}

}  // namespace fluoroseq
