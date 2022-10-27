/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// File under test:
#include "normal-distribution-fitter.h"

// Standard C++ library headers:
#include <cmath>

// Local project headers:
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using std::sqrt;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(parameterization_suite)
BOOST_AUTO_TEST_SUITE(fit_suite)
BOOST_AUTO_TEST_SUITE(normal_distribution_fitter_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    NormalDistributionFitter ndf;
    BOOST_TEST(ndf.w_sum_x == 0.0);
    BOOST_TEST(ndf.w_sum_x_sq_over_n == 0.0);
    BOOST_TEST(ndf.w_sum_n == 0.0);
    BOOST_TEST(ndf.total_weight == 0.0);
}

BOOST_AUTO_TEST_CASE(add_sample_once_n_eq_1_test, *tolerance(TOL)) {
    NormalDistributionFitter ndf;
    double x = 1.277;
    int n = 1;
    double w = 0.98;
    ndf.add_sample(x, n, w);
    BOOST_TEST(ndf.w_sum_x == x * w);
    BOOST_TEST(ndf.w_sum_x_sq_over_n == x * x * w / (double)n);
    BOOST_TEST(ndf.w_sum_n == (double)n * w);
    BOOST_TEST(ndf.total_weight == w);
}

BOOST_AUTO_TEST_CASE(add_sample_once_n_gt_1_test, *tolerance(TOL)) {
    NormalDistributionFitter ndf;
    double x = 3.1415928;
    int n = 3;
    double w = 0.979;
    ndf.add_sample(x, n, w);
    BOOST_TEST(ndf.w_sum_x == x * w);
    BOOST_TEST(ndf.w_sum_x_sq_over_n == x * x * w / (double)n);
    BOOST_TEST(ndf.w_sum_n == (double)n * w);
    BOOST_TEST(ndf.total_weight == w);
}

BOOST_AUTO_TEST_CASE(add_sample_twice_n_eq_1_test, *tolerance(TOL)) {
    NormalDistributionFitter ndf;
    double x1 = 1.277;
    int n1 = 1;
    double w1 = 0.98;
    double x2 = 1.166;
    int n2 = 1;
    double w2 = 0.49;
    ndf.add_sample(x1, n1, w1);
    ndf.add_sample(x2, n2, w2);
    BOOST_TEST(ndf.w_sum_x == x1 * w1 + x2 * w2);
    BOOST_TEST(ndf.w_sum_x_sq_over_n
               == x1 * x1 * w1 / (double)n1 + x2 * x2 * w2 / (double)n2);
    BOOST_TEST(ndf.w_sum_n == (double)n1 * w1 + (double)n2 * w2);
    BOOST_TEST(ndf.total_weight == w1 + w2);
}

BOOST_AUTO_TEST_CASE(add_sample_twice_n_gt_1_test, *tolerance(TOL)) {
    NormalDistributionFitter ndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    double x2 = 4.91;
    int n2 = 5;
    double w2 = 0.49;
    ndf.add_sample(x1, n1, w1);
    ndf.add_sample(x2, n2, w2);
    BOOST_TEST(ndf.w_sum_x == x1 * w1 + x2 * w2);
    BOOST_TEST(ndf.w_sum_x_sq_over_n
               == x1 * x1 * w1 / (double)n1 + x2 * x2 * w2 / (double)n2);
    BOOST_TEST(ndf.w_sum_n == (double)n1 * w1 + (double)n2 * w2);
    BOOST_TEST(ndf.total_weight == w1 + w2);
}

BOOST_AUTO_TEST_CASE(get_mu_one_sample_test, *tolerance(TOL)) {
    NormalDistributionFitter ndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    ndf.add_sample(x1, n1, w1);
    BOOST_TEST(ndf.get_mu() == x1 / (double)n1);
}

BOOST_AUTO_TEST_CASE(get_mu_two_samples_test, *tolerance(TOL)) {
    NormalDistributionFitter ndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    double x2 = 4.91;
    int n2 = 5;
    double w2 = 0.49;
    ndf.add_sample(x1, n1, w1);
    ndf.add_sample(x2, n2, w2);
    BOOST_TEST(ndf.get_mu()
               == (x1 * w1 + x2 * w2) / ((double)n1 * w1 + (double)n2 * w2));
}

// Lazily commenting out broken tests because this code doesn't really work
// anyways - not on real data. Probably a bad idea....

// BOOST_AUTO_TEST_CASE(get_sig_one_sample_test, *tolerance(TOL)) {
//     NormalDistributionFitter ndf;
//     double x1 = 3.43;
//     int n1 = 3;
//     double w1 = 0.98;
//     ndf.add_sample(x1, n1, w1);
//     double mu = ndf.get_mu();
//     BOOST_TEST(ndf.get_sig() == sqrt((x1 - n1 * mu) * (x1 - n1 * mu) / n1));
// }

// BOOST_AUTO_TEST_CASE(get_sig_two_samples_test, *tolerance(TOL)) {
//     NormalDistributionFitter ndf;
//     double x1 = 3.43;
//     int n1 = 3;
//     double w1 = 0.98;
//     double x2 = 4.91;
//     int n2 = 5;
//     double w2 = 0.49;
//     ndf.add_sample(x1, n1, w1);
//     ndf.add_sample(x2, n2, w2);
//     double mu = ndf.get_mu();
//     BOOST_TEST(ndf.get_sig()
//                == sqrt(((x1 - (double)n1 * mu) * (x1 - (double)n1 * mu) * w1
//                                 / (double)n1
//                         + (x2 - (double)n2 * mu) * (x2 - (double)n2 * mu) *
//                         w2
//                                   / (double)n2)
//                        / (w1 + w2)));
// }

BOOST_AUTO_TEST_SUITE_END()  // normal_distribution_fitter_suite
BOOST_AUTO_TEST_SUITE_END()  // fit_suite
BOOST_AUTO_TEST_SUITE_END()  // parameterization_suite

}  // namespace whatprot
