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
#include "error-model.h"

// Standard C++ library headers:
#include <cmath>
#include <functional>

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using std::exp;
using std::function;
using std::log;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(common_suite)
BOOST_AUTO_TEST_SUITE(error_model_suite)

BOOST_AUTO_TEST_CASE(constructor_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    BOOST_TEST(em.p_edman_failure == p_edman_failure);
    BOOST_TEST(em.p_detach == p_detach);
    BOOST_TEST(em.p_bleach == p_bleach);
    BOOST_TEST(em.p_dud == p_dud);
    BOOST_TEST(em.distribution_type == dist_type);
    BOOST_TEST(em.mu == mu);
    BOOST_TEST(em.sigma == sigma);
    BOOST_TEST(em.stuck_dye_ratio == stuck_dye_ratio);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_zero_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 0.0;
    int state = 0;
    BOOST_TEST(pdf(observed, state) == 1.0);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_zero_obs_one_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 1.0;
    int state = 0;
    BOOST_TEST(pdf(observed, state) == 0.0);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_one_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 0.0;
    int state = 1;
    BOOST_TEST(pdf(observed, state) == 0.0);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_one_obs_one_test, *tolerance(TOL)) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 1.0;
    int state = 1;
    // The test value was found using an online lognormal distribution pdf
    // calculator.
    BOOST_TEST(pdf(observed, state) == 2.4933892525089547);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_eq_obs_ne_one_test, *tolerance(TOL)) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.3);
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 1.3;
    int state = 1;
    // The test value was found using an online lognormal distribution pdf
    // calculator.
    BOOST_TEST(pdf(observed, state) == 2.4933892525089547 / 1.3);
}

BOOST_AUTO_TEST_CASE(pdf_override_state_zero_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 0.0;
    int state = 0;
    BOOST_TEST(pdf(observed, state) == 1.0);
}

BOOST_AUTO_TEST_CASE(pdf_override_state_zero_obs_one_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 1.0;
    int state = 0;
    BOOST_TEST(pdf(observed, state) == 1.0);
}

BOOST_AUTO_TEST_CASE(pdf_override_state_one_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 0.0;
    int state = 1;
    BOOST_TEST(pdf(observed, state) == 1.0);
}

BOOST_AUTO_TEST_CASE(pdf_override_state_one_obs_one_test, *tolerance(TOL)) {
    double p_edman_failure = .07;
    double p_detach = .04;
    double p_bleach = .05;
    double p_dud = .10;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = .16;
    double stuck_dye_ratio = 0.5;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio);
    function<double(double, int)> pdf = em.pdf();
    double observed = 1.0;
    int state = 1;
    BOOST_TEST(pdf(observed, state) == 1.0);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_edman_failure_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    ErrorModel em2(
            0.66, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    BOOST_TEST(em1.relative_distance(em2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(em2.relative_distance(em1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_detach_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    ErrorModel em2(
            0.5, 0.66, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    BOOST_TEST(em1.relative_distance(em2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(em2.relative_distance(em1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_bleach_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    ErrorModel em2(
            0.5, 0.5, 0.66, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    BOOST_TEST(em1.relative_distance(em2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(em2.relative_distance(em1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_dud_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    ErrorModel em2(
            0.5, 0.5, 0.5, 0.66, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    BOOST_TEST(em1.relative_distance(em2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(em2.relative_distance(em1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_mu_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    ErrorModel em2(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.66, 0.5, 0.5);
    BOOST_TEST(em1.relative_distance(em2) == (exp(0.66) - exp(0.5)) / exp(0.5));
    BOOST_TEST(em2.relative_distance(em1)
               == (exp(0.66) - exp(0.5)) / exp(0.66));
}

BOOST_AUTO_TEST_CASE(relative_distance_sigma_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    ErrorModel em2(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.66, 0.5);
    BOOST_TEST(em1.relative_distance(em2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(em2.relative_distance(em1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_stuck_dye_ratio_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.5);
    ErrorModel em2(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.5, 0.5, 0.66);
    BOOST_TEST(em1.relative_distance(em2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(em2.relative_distance(em1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_max_no_sum_test, *tolerance(TOL)) {
    ErrorModel em1(
            0.5, 0.5, 0.5, 0.5, DistributionType::OVERRIDE, 0.66, 0.5, 0.5);
    ErrorModel em2(
            0.7, 0.7, 0.7, 0.7, DistributionType::OVERRIDE, 0.66, 0.7, 0.7);
    BOOST_TEST(em1.relative_distance(em2) == (0.7 - 0.5) / 0.5);
    BOOST_TEST(em2.relative_distance(em1) == (0.7 - 0.5) / 0.7);
}

BOOST_AUTO_TEST_SUITE_END()  // error_model_suite
BOOST_AUTO_TEST_SUITE_END()  // common_suite

}  // namespace whatprot
