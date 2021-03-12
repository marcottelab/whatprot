/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "normal-distribution-fitter.h"

// Standard C++ library headers:
#include <cmath>

// Local project headers:
#include "common/error-model.h"

namespace whatprot {

namespace {
using std::sqrt;
}  // namespace

NormalDistributionFitter::NormalDistributionFitter()
        : w_sum_x(0.0),
          w_sum_x_sq_over_n(0.0),
          w_sum_n(0.0),
          total_weight(0.0) {}

void NormalDistributionFitter::add_sample(double x, int n, double weight) {
    if (n == 0) {
        return;
    }
    w_sum_x += x * weight;
    // For estimating the variance, x needs to be scaled by the square root of
    // n. This can be verified analytically. Note that this was not needed
    // for the w_sum_x variable (again, can be verified analytically).
    w_sum_x_sq_over_n += x * (x * weight) / (double)n;
    w_sum_n += (double)n * weight;
    total_weight += weight;
}

DistributionType NormalDistributionFitter::get_type() const {
    return DistributionType::NORMAL;
}

double NormalDistributionFitter::get_mu() const {
    return w_sum_x / w_sum_n;
}

double NormalDistributionFitter::get_sigma() const {
    double mu = get_mu();
    double sigma_sq = (w_sum_x_sq_over_n - mu * mu * w_sum_n) / total_weight;
    return sqrt(sigma_sq);
}

NormalDistributionFitter NormalDistributionFitter::operator+(
        const NormalDistributionFitter& other) const {
    NormalDistributionFitter result_fitter;
    result_fitter.w_sum_x = w_sum_x + other.w_sum_x;
    result_fitter.w_sum_x_sq_over_n = w_sum_x_sq_over_n + other.w_sum_x_sq_over_n;
    result_fitter.w_sum_n = w_sum_n + other.w_sum_n;
    result_fitter.total_weight = total_weight + other.total_weight;
    return result_fitter;
}

void NormalDistributionFitter::operator+=(
        const NormalDistributionFitter& other) {
    w_sum_x += other.w_sum_x;
    w_sum_x_sq_over_n += other.w_sum_x_sq_over_n;
    w_sum_n += other.w_sum_n;
    total_weight += other.total_weight;
}

void NormalDistributionFitter::operator*=(double weight_adjustment) {
    w_sum_x *= weight_adjustment;
    w_sum_x_sq_over_n *= weight_adjustment;
    w_sum_n *= weight_adjustment;
    total_weight *= weight_adjustment;
}

}  // namespace whatprot
