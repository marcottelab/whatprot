// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_COMMON_RADIOMETRY_H
#define FLUOROSEQ_COMMON_RADIOMETRY_H

namespace fluoroseq {

class Radiometry {
public:
    Radiometry(int num_timesteps, int num_channels);
    Radiometry(const Radiometry& other);
    ~Radiometry();
    double& operator()(int t, int c);
    double operator()(int t, int c) const;

    double* intensities;
    int num_timesteps;
    int num_channels;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_RADIOMETRY_H
