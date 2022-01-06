/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_FIT_LOG_NORMAL_DISTRIBUTION_FITTER_H
#define WHATPROT_PARAMETERIZATION_FIT_LOG_NORMAL_DISTRIBUTION_FITTER_H

// Local project headers:
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

class LogNormalDistributionFitter {
public:
    LogNormalDistributionFitter();
    virtual ~LogNormalDistributionFitter();
    virtual void add_sample(double x, int n, double weight);
    double get_mu() const;
    double get_sig() const;
    LogNormalDistributionFitter operator+(
            const LogNormalDistributionFitter& other) const;
    void operator+=(const LogNormalDistributionFitter& other);
    void operator*=(double weight_adjustment);
    double w_sum_log_x_over_n;
    double w_sum_log_x_over_n_sq;
    double total_weight;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_FIT_LOG_NORMAL_DISTRIBUTION_FITTER_H
