/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_COMMON_RADIOMETRY_H
#define WHATPROT_COMMON_RADIOMETRY_H

namespace whatprot {

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

}  // namespace whatprot

#endif  // WHATPROT_COMMON_RADIOMETRY_H
