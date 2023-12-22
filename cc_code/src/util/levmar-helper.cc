/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "levmar-helper.h"

// Standard C++ library headers:
#include <cmath>
#include <vector>

// External headers:
#include "levmar.h"

namespace {
using std::exp;
using std::sqrt;
using std::vector;
}  // namespace

namespace whatprot {

void offset_exponential(double* p, double* y, int m, int n, void* data) {
    double* sqrt_w = (double*)data;
    for (int i = 0; i < n; i++) {
        y[i] = sqrt_w[i] * (p[0] + p[1] * exp(-p[2] * i));
    }
}

void jacobian_of_offset_exponential(
        double* p, double* jac, int m, int n, void* data) {
    double* sqrt_w = (double*)data;
    unsigned int j = 0;
    for (int i = 0; i < n; i++) {
        double tmp = exp(-p[2] * i);
        jac[j++] = sqrt_w[i];
        jac[j++] = sqrt_w[i] * tmp;
        jac[j++] = sqrt_w[i] * -p[1] * i * tmp;
    }
}

void least_squares_fit_of_offset_exponential(const vector<double>& y,
                                             const vector<double>& w,
                                             const vector<bool>& hold,
                                             vector<double>* p) {
    // We use bounds to hold parameters constant.
    vector<double> lb(3, 0);  // lower bound.
    vector<double> ub(3, 1);  // upper bound.
    for (unsigned int i = 0; i < 3; i++) {
        if (hold[i] == true) {
            lb[i] = (*p)[i];
            ub[i] = (*p)[i];
        }
    }
    // Least squares fit with weights changes X'XA=X'Y to X'WXA=X'WY, where W
    // is the weight matrix, with weights along the diagonal. This
    // implementation of least squares fit does not accept a weight matrix, but
    // by multiplying X (during function evaluation) and Y by sqrt(W) we get
    // the same behavior.
    vector<double> sqrt_w(w.size());
    for (unsigned int i = 0; i < w.size(); i++) {
        sqrt_w[i] = sqrt(w[i]);
    }
    vector<double> weighted_y(y.size());
    for (unsigned int i = 0; i < y.size(); i++) {
        weighted_y[i] = sqrt_w[i] * y[i];
    }
    dlevmar_bc_der(&offset_exponential,  // func
                   &jacobian_of_offset_exponential,  // jacf
                   &(*p)[0],  // p
                   &weighted_y[0],  // "x" (official name - y is a better name)
                   p->size(),  // m
                   weighted_y.size(),  // n
                   &lb[0],  // lb
                   &ub[0],  // ub
                   NULL,  // dscl
                   10,  // itmax
                   NULL,  // opts
                   NULL,  // info
                   NULL,  // work
                   NULL,  // covar
                   &sqrt_w[0]);  // adata (any extra data we want to send)
}

}  // namespace whatprot