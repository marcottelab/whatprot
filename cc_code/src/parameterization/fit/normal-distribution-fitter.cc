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
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using std::sqrt;
}  // namespace

NormalDistributionFitter::NormalDistributionFitter()
        : w_sum_x(0.0),
          w_sum_x_sq_over_n(0.0),
          w_sum_n(0.0),
          total_weight(0.0) {}

NormalDistributionFitter::~NormalDistributionFitter() {}

void NormalDistributionFitter::add_sample(double x, int n, double weight) {
    if (n == 0) {
        zero_w_sum_x_sq += x * (x * weight);
        zero_total_weight += weight;
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

double NormalDistributionFitter::get_mu() const {
    return w_sum_x / w_sum_n;
}

double NormalDistributionFitter::get_bg_sig() const {
    // double bg_sig_sq = zero_w_sum_x_sq / zero_total_weight;
    // return sqrt(bg_sig_sq);
    return .07;
}

double NormalDistributionFitter::get_sig() const {
    double mu = get_mu();
    double bg_sig = get_bg_sig();
    double sig_sq = (w_sum_x_sq_over_n - mu * mu * w_sum_n) / total_weight;
    return sqrt(sig_sq - bg_sig * bg_sig);
}

NormalDistributionFitter NormalDistributionFitter::operator+(
        const NormalDistributionFitter& other) const {
    NormalDistributionFitter result_fitter;
    result_fitter.w_sum_x = w_sum_x + other.w_sum_x;
    result_fitter.w_sum_x_sq_over_n =
            w_sum_x_sq_over_n + other.w_sum_x_sq_over_n;
    result_fitter.w_sum_n = w_sum_n + other.w_sum_n;
    result_fitter.total_weight = total_weight + other.total_weight;
    result_fitter.zero_w_sum_x_sq = zero_w_sum_x_sq + other.zero_w_sum_x_sq;
    result_fitter.zero_total_weight =
            zero_total_weight + other.zero_total_weight;
    return result_fitter;
}

void NormalDistributionFitter::operator+=(
        const NormalDistributionFitter& other) {
    w_sum_x += other.w_sum_x;
    w_sum_x_sq_over_n += other.w_sum_x_sq_over_n;
    w_sum_n += other.w_sum_n;
    total_weight += other.total_weight;
    zero_w_sum_x_sq += other.zero_w_sum_x_sq;
    zero_total_weight += other.zero_total_weight;
}

void NormalDistributionFitter::operator*=(double weight_adjustment) {
    w_sum_x *= weight_adjustment;
    w_sum_x_sq_over_n *= weight_adjustment;
    w_sum_n *= weight_adjustment;
    total_weight *= weight_adjustment;
    zero_w_sum_x_sq *= weight_adjustment;
    zero_w_sum_x_sq *= weight_adjustment;
}

}  // namespace whatprot
