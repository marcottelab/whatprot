/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_FIT_NORMAL_DISTRIBUTION_FITTER_H
#define WHATPROT_PARAMETERIZATION_FIT_NORMAL_DISTRIBUTION_FITTER_H

// Local project headers:
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

class NormalDistributionFitter {
public:
    NormalDistributionFitter();
    virtual ~NormalDistributionFitter();
    virtual void add_sample(double x, int n, double weight);
    double get_mu() const;
    double get_sigma() const;
    NormalDistributionFitter operator+(
            const NormalDistributionFitter& other) const;
    void operator+=(const NormalDistributionFitter& other);
    void operator*=(double weight_adjustment);
    double w_sum_x;
    double w_sum_x_sq_over_n;
    double w_sum_n;
    double total_weight;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_FIT_NORMAL_DISTRIBUTION_FITTER_H
