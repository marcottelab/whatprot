/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_FIT_LOG_NORMAL_DISTRIBUTION_FITTER_H
#define WHATPROT_HMM_FIT_LOG_NORMAL_DISTRIBUTION_FITTER_H

// Local project headers:
#include "common/error-model.h"

namespace whatprot {

class LogNormalDistributionFitter {
public:
    LogNormalDistributionFitter();
    void add_sample(double x, int n, double weight);
    DistributionType get_type() const;
    double get_mu() const;
    double get_sigma() const;
    LogNormalDistributionFitter operator+(
            const LogNormalDistributionFitter& other) const;
    void operator+=(const LogNormalDistributionFitter& other);
    void operator*=(double weight_adjustment);
    double w_sum_log_x_over_n;
    double w_sum_log_x_over_n_sq;
    double total_weight;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_FIT_LOG_NORMAL_DISTRIBUTION_FITTER_H
