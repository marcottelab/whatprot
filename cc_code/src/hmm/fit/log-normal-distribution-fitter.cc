/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "log-normal-distribution-fitter.h"

// Standard C++ library headers:
#include <cmath>

// Local project headers:
#include "common/error-model.h"

namespace whatprot {

namespace {
using std::log;
using std::sqrt;
}  // namespace

LogNormalDistributionFitter::LogNormalDistributionFitter()
        : w_sum_log_x_over_n(0.0),
          w_sum_log_x_over_n_sq(0.0),
          total_weight(0.0) {}

void LogNormalDistributionFitter::add_sample(double x, int n, double weight) {
    if (n == 0 || x == 0.0) {
        return;
    }
    // Need to divide x by n before taking log to get the scaling right. This
    // scaling can be verified analytically. Once we've done this the Maximum
    // Likelihood Estimator can be found simply by finding the normal
    // distribution that best fits these values.
    double log_x_over_n = log(x / (double)n);
    w_sum_log_x_over_n += weight * log_x_over_n;
    w_sum_log_x_over_n_sq += weight * log_x_over_n * log_x_over_n;
    total_weight += weight;
}

DistributionType LogNormalDistributionFitter::get_type() const {
    return DistributionType::LOGNORMAL;
}

double LogNormalDistributionFitter::get_mu() const {
    return w_sum_log_x_over_n / total_weight;
}

double LogNormalDistributionFitter::get_sigma() const {
    double mu = get_mu();
    return sqrt(w_sum_log_x_over_n_sq / total_weight - mu * mu);
}

LogNormalDistributionFitter LogNormalDistributionFitter::operator+(
        const LogNormalDistributionFitter& other) const {
    LogNormalDistributionFitter result_fitter;
    result_fitter.w_sum_log_x_over_n =
            w_sum_log_x_over_n + other.w_sum_log_x_over_n;
    result_fitter.w_sum_log_x_over_n_sq =
            w_sum_log_x_over_n_sq + other.w_sum_log_x_over_n_sq;
    result_fitter.total_weight = total_weight + other.total_weight;
    return result_fitter;
}

void LogNormalDistributionFitter::operator+=(
        const LogNormalDistributionFitter& other) {
    w_sum_log_x_over_n += other.w_sum_log_x_over_n;
    w_sum_log_x_over_n_sq += other.w_sum_log_x_over_n_sq;
    total_weight += other.total_weight;
}

void LogNormalDistributionFitter::operator*=(double weight_adjustment) {
    w_sum_log_x_over_n *= weight_adjustment;
    w_sum_log_x_over_n_sq *= weight_adjustment;
    total_weight *= weight_adjustment;
}

}  // namespace whatprot
