/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "error-model.h"

// Standard C++ library headers:
#include <algorithm>
#include <cmath>
#include <functional>
#include <string>

namespace whatprot {

namespace {
using std::abs;
using std::exp;
using std::function;
using std::log;
using std::max;
using std::sqrt;
using std::string;
using std::to_string;
double PI = 3.141592653589793238;
}  // namespace

ErrorModel::ErrorModel(double p_edman_failure,
                       double p_detach,
                       double p_bleach,
                       double p_dud,
                       DistributionType distribution_type,
                       double mu,
                       double sigma,
                       double stuck_dye_ratio)
        : p_edman_failure(p_edman_failure),
          p_detach(p_detach),
          p_bleach(p_bleach),
          p_dud(p_dud),
          distribution_type(distribution_type),
          mu(mu),
          sigma(sigma),
          stuck_dye_ratio(stuck_dye_ratio) {}

function<double(double, int)> ErrorModel::pdf() const {
    switch (distribution_type) {
        case OVERRIDE:
            return [](double observed, int state) -> double {
                return 1.0;
            };
            break;
        case LOGNORMAL:
        default:
            double mu = this->mu;
            double sigma = this->sigma;
            double multiplier = 1.0 / (sigma * sqrt(2.0 * PI));
            return [mu, sigma, multiplier](double observed,
                                           int state) -> double {
                if (state > 0) {
                    if (observed == 0.0) {
                        return 0.0;
                    } else {
                        double offset = log(observed) - log((double)state) - mu;
                        return (multiplier / observed)
                               * exp(-(offset * offset)
                                     / (2.0 * sigma * sigma));
                    }
                } else {
                    if (observed == 0.0) {
                        return 1.0;
                    } else {
                        return 0.0;
                    }
                }
            };
    }
}

double ErrorModel::relative_distance(const ErrorModel& error_model) const {
    double dist = 0.0;
    dist = max(dist,
               abs(p_edman_failure - error_model.p_edman_failure)
                       / p_edman_failure);
    dist = max(dist, abs(p_detach - error_model.p_detach) / p_detach);
    dist = max(dist, abs(p_bleach - error_model.p_bleach) / p_bleach);
    dist = max(dist, abs(p_dud - error_model.p_dud) / p_dud);
    dist = max(dist, abs(exp(mu) - exp(error_model.mu)) / exp(mu));
    dist = max(dist, abs(sigma - error_model.sigma) / sigma);
    dist = max(dist,
               abs(stuck_dye_ratio - error_model.stuck_dye_ratio)
                       / stuck_dye_ratio);
    return dist;
}

string ErrorModel::debug_string() const {
    return "Edman failure rate: " + to_string(p_edman_failure)
           + ", Detach rate: " + to_string(p_detach) + ", Bleach rate: "
           + to_string(p_bleach) + ", Dud rate: " + to_string(p_dud)
           + ", exp(mu): " + to_string(exp(mu)) + ", sigma: " + to_string(sigma)
           + ", Stuck dye ratio: " + to_string(stuck_dye_ratio);
}

}  // namespace whatprot
