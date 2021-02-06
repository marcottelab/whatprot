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
#include "log-normal-distribution-fitter.h"

// Standard C++ library headers:
#include <cmath>

// Local project headers:
#include "common/error-model.h"

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
using std::log;
using std::sqrt;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(fit_suite)
BOOST_AUTO_TEST_SUITE(log_normal_distribution_fitter_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    BOOST_TEST(lndf.w_sum_log_x_over_n == 0.0);
    BOOST_TEST(lndf.w_sum_log_x_over_n_sq == 0.0);
    BOOST_TEST(lndf.total_weight == 0.0);
}

BOOST_AUTO_TEST_CASE(add_sample_once_n_eq_1_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x = 1.277;
    int n = 1;
    double w = 0.98;
    lndf.add_sample(x, n, w);
    BOOST_TEST(lndf.w_sum_log_x_over_n == log(x / (double)n) * w);
    BOOST_TEST(lndf.w_sum_log_x_over_n_sq
               == log(x / (double)n) * log(x / (double)n) * w);
    BOOST_TEST(lndf.total_weight == w);
}

BOOST_AUTO_TEST_CASE(add_sample_once_n_gt_1_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x = 3.1415928;
    int n = 3;
    double w = 0.979;
    lndf.add_sample(x, n, w);
    BOOST_TEST(lndf.w_sum_log_x_over_n == log(x / (double)n) * w);
    BOOST_TEST(lndf.w_sum_log_x_over_n_sq
               == log(x / (double)n) * log(x / (double)n) * w);
    BOOST_TEST(lndf.total_weight == w);
}

BOOST_AUTO_TEST_CASE(add_sample_twice_n_eq_1_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x1 = 1.277;
    int n1 = 1;
    double w1 = 0.98;
    double x2 = 1.166;
    int n2 = 1;
    double w2 = 0.49;
    lndf.add_sample(x1, n1, w1);
    lndf.add_sample(x2, n2, w2);
    BOOST_TEST(lndf.w_sum_log_x_over_n
               == log(x1 / (double)n1) * w1 + log(x2 / (double)n2) * w2);
    BOOST_TEST(lndf.w_sum_log_x_over_n_sq
               == log(x1 / (double)n1) * log(x1 / (double)n1) * w1
                          + log(x2 / (double)n2) * log(x2 / (double)n2) * w2);
    BOOST_TEST(lndf.total_weight == w1 + w2);
}

BOOST_AUTO_TEST_CASE(add_sample_twice_n_gt_1_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    double x2 = 4.91;
    int n2 = 5;
    double w2 = 0.49;
    lndf.add_sample(x1, n1, w1);
    lndf.add_sample(x2, n2, w2);
    BOOST_TEST(lndf.w_sum_log_x_over_n
               == log(x1 / (double)n1) * w1 + log(x2 / (double)n2) * w2);
    BOOST_TEST(lndf.w_sum_log_x_over_n_sq
               == log(x1 / (double)n1) * log(x1 / (double)n1) * w1
                          + log(x2 / (double)n2) * log(x2 / (double)n2) * w2);
    BOOST_TEST(lndf.total_weight == w1 + w2);
}

BOOST_AUTO_TEST_CASE(get_type_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    BOOST_TEST(lndf.get_type() == DistributionType::LOGNORMAL);
}

BOOST_AUTO_TEST_CASE(get_mu_one_sample_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    lndf.add_sample(x1, n1, w1);
    BOOST_TEST(lndf.get_mu() == log(x1 / (double)n1));
}

BOOST_AUTO_TEST_CASE(get_mu_two_samples_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    double x2 = 4.91;
    int n2 = 5;
    double w2 = 0.49;
    lndf.add_sample(x1, n1, w1);
    lndf.add_sample(x2, n2, w2);
    BOOST_TEST(lndf.get_mu()
               == (log(x1 / (double)n1) * w1 + log(x2 / (double)n2) * w2)
                          / (w1 + w2));
}

BOOST_AUTO_TEST_CASE(get_sigma_one_sample_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    lndf.add_sample(x1, n1, w1);
    double mu = lndf.get_mu();
    BOOST_TEST(lndf.get_sigma()
               == sqrt((log(x1 / n1) - mu) * (log(x1 / n1) - mu)));
}

BOOST_AUTO_TEST_CASE(get_sigma_two_samples_test, *tolerance(TOL)) {
    LogNormalDistributionFitter lndf;
    double x1 = 3.43;
    int n1 = 3;
    double w1 = 0.98;
    double x2 = 4.91;
    int n2 = 5;
    double w2 = 0.49;
    lndf.add_sample(x1, n1, w1);
    lndf.add_sample(x2, n2, w2);
    double mu = lndf.get_mu();
    BOOST_TEST(lndf.get_sigma()
               == sqrt(((log(x1 / n1) - mu) * (log(x1 / n1) - mu) * w1
                        + (log(x2 / n2) - mu) * (log(x2 / n2) - mu) * w2)
                       / (w1 + w2)));
}

BOOST_AUTO_TEST_SUITE_END()  // log_normal_distribution_fitter_suite
BOOST_AUTO_TEST_SUITE_END()  // fit_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace fluoroseq
