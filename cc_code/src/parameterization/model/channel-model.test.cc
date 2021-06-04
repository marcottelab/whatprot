/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin *
* Department: Oden Institute and Institute for Cellular and Molecular Biology
*
* PI: Edward Marcotte *
* Project: Protein Fluorosequencing *
\******************************************************************************/

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// File under test:
#include "channel-model.h"

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

BOOST_AUTO_TEST_SUITE(parameterization_suite)
BOOST_AUTO_TEST_SUITE(model_suite)
BOOST_AUTO_TEST_SUITE(channel_model_suite)

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_zero_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.mu = log(1.0);
    channel_model.sigma = .16;
    channel_model.stuck_dye_ratio = 0.5;
    channel_model.p_stuck_dye_loss = 0.08;
    double observed = 0.0;
    int state = 0;
    BOOST_TEST(channel_model.pdf(observed, state) == 1.0);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_zero_obs_one_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.mu = log(1.0);
    channel_model.sigma = .16;
    channel_model.stuck_dye_ratio = 0.5;
    channel_model.p_stuck_dye_loss = 0.08;
    double observed = 1.0;
    int state = 0;
    BOOST_TEST(channel_model.pdf(observed, state) == 0.0);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_one_obs_zero_test) {
    double p_edman_failure = .07;
    double p_detach = .04;
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.mu = log(1.0);
    channel_model.sigma = .16;
    channel_model.stuck_dye_ratio = 0.5;
    channel_model.p_stuck_dye_loss = 0.08;
    double observed = 0.0;
    int state = 1;
    BOOST_TEST(channel_model.pdf(observed, state) == 0.0);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_one_obs_one_test, *tolerance(TOL)) {
    double p_edman_failure = .07;
    double p_detach = .04;
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.mu = log(1.0);
    channel_model.sigma = .16;
    channel_model.stuck_dye_ratio = 0.5;
    channel_model.p_stuck_dye_loss = 0.08;
    double observed = 1.0;
    int state = 1;
    // The test value was found using an online lognormal distribution pdf
    // calculator.
    BOOST_TEST(channel_model.pdf(observed, state) == 2.4933892525089547);
}

BOOST_AUTO_TEST_CASE(pdf_lognormal_state_eq_obs_ne_one_test, *tolerance(TOL)) {
    double p_edman_failure = .07;
    double p_detach = .04;
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.mu = log(1.3);
    channel_model.sigma = .16;
    channel_model.stuck_dye_ratio = 0.5;
    channel_model.p_stuck_dye_loss = 0.08;
    double observed = 1.3;
    int state = 1;
    // The test value was found using an online lognormal distribution pdf
    // calculator.
    BOOST_TEST(channel_model.pdf(observed, state) == 2.4933892525089547 / 1.3);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_bleach_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.mu = 0.5;
    cm1.sigma = 0.5;
    cm1.stuck_dye_ratio = 0.5;
    cm1.p_stuck_dye_loss = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.66;
    cm2.p_dud = 0.5;
    cm2.mu = 0.5;
    cm2.sigma = 0.5;
    cm2.stuck_dye_ratio = 0.5;
    cm2.p_stuck_dye_loss = 0.5;

    BOOST_TEST(cm1.relative_distance(cm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(cm2.relative_distance(cm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_dud_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.mu = 0.5;
    cm1.sigma = 0.5;
    cm1.stuck_dye_ratio = 0.5;
    cm1.p_stuck_dye_loss = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.66;
    cm2.mu = 0.5;
    cm2.sigma = 0.5;
    cm2.stuck_dye_ratio = 0.5;
    cm2.p_stuck_dye_loss = 0.5;

    BOOST_TEST(cm1.relative_distance(cm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(cm2.relative_distance(cm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_mu_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.mu = 0.5;
    cm1.sigma = 0.5;
    cm1.stuck_dye_ratio = 0.5;
    cm1.p_stuck_dye_loss = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.5;
    cm2.mu = 0.66;
    cm2.sigma = 0.5;
    cm2.stuck_dye_ratio = 0.5;
    cm2.p_stuck_dye_loss = 0.5;

    BOOST_TEST(cm1.relative_distance(cm2) == (exp(0.66) - exp(0.5)) / exp(0.5));
    BOOST_TEST(cm2.relative_distance(cm1)
               == (exp(0.66) - exp(0.5)) / exp(0.66));
}

BOOST_AUTO_TEST_CASE(relative_distance_sigma_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.mu = 0.5;
    cm1.sigma = 0.5;
    cm1.stuck_dye_ratio = 0.5;
    cm1.p_stuck_dye_loss = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.5;
    cm2.mu = 0.5;
    cm2.sigma = 0.66;
    cm2.stuck_dye_ratio = 0.5;
    cm2.p_stuck_dye_loss = 0.5;

    BOOST_TEST(cm1.relative_distance(cm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(cm2.relative_distance(cm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_stuck_dye_ratio_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.mu = 0.5;
    cm1.sigma = 0.5;
    cm1.stuck_dye_ratio = 0.5;
    cm1.p_stuck_dye_loss = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.5;
    cm2.mu = 0.5;
    cm2.sigma = 0.5;
    cm2.stuck_dye_ratio = 0.66;
    cm2.p_stuck_dye_loss = 0.5;

    BOOST_TEST(cm1.relative_distance(cm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(cm2.relative_distance(cm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_stuck_dye_loss_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.mu = 0.5;
    cm1.sigma = 0.5;
    cm1.stuck_dye_ratio = 0.5;
    cm1.p_stuck_dye_loss = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.5;
    cm2.mu = 0.5;
    cm2.sigma = 0.5;
    cm2.stuck_dye_ratio = 0.5;
    cm2.p_stuck_dye_loss = 0.66;

    BOOST_TEST(cm1.relative_distance(cm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(cm2.relative_distance(cm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_max_no_sum_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.mu = 0.66;
    cm1.sigma = 0.5;
    cm1.stuck_dye_ratio = 0.5;
    cm1.p_stuck_dye_loss = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.7;
    cm2.p_dud = 0.7;
    cm2.mu = 0.66;
    cm2.sigma = 0.7;
    cm2.stuck_dye_ratio = 0.7;
    cm2.p_stuck_dye_loss = 0.7;

    BOOST_TEST(cm1.relative_distance(cm2) == (0.7 - 0.5) / 0.5);
    BOOST_TEST(cm2.relative_distance(cm1) == (0.7 - 0.5) / 0.7);
}

BOOST_AUTO_TEST_SUITE_END()  // channel_model_suite
BOOST_AUTO_TEST_SUITE_END()  // model_suite
BOOST_AUTO_TEST_SUITE_END()  // parameterization_suite

}  // namespace whatprot
