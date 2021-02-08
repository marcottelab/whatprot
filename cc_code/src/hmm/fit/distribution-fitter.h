/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_FIT_DISTRIBUTION_FITTER_H
#define WHATPROT_HMM_FIT_DISTRIBUTION_FITTER_H

// Local project headers:
#include "common/error-model.h"

namespace whatprot {

class DistributionFitter {
public:
    virtual void add_sample(double x, int n, double weight) = 0;
    virtual DistributionType get_type() const = 0;
    virtual double get_mu() const = 0;
    virtual double get_sigma() const = 0;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_FIT_DISTRIBUTION_FITTER_H
