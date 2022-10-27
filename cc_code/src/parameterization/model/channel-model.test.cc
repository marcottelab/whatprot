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
using std::sqrt;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(parameterization_suite)
BOOST_AUTO_TEST_SUITE(model_suite)
BOOST_AUTO_TEST_SUITE(channel_model_suite)

BOOST_AUTO_TEST_CASE(pdf_state_zero_obs_zero_test) {
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 0.0;
    int state = 0;
    BOOST_TEST(channel_model.pdf(observed, state) == 59.811436342043883);
}

BOOST_AUTO_TEST_CASE(pdf_state_zero_obs_one_test) {
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 1.0;
    int state = 0;
    BOOST_TEST(channel_model.pdf(observed, state) == 0);
}

BOOST_AUTO_TEST_CASE(pdf_state_one_obs_zero_test) {
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 0.0;
    int state = 1;
    BOOST_TEST(channel_model.pdf(observed, state) == 8.488175272749065e-09);
}

BOOST_AUTO_TEST_CASE(pdf_state_one_obs_one_test, *tolerance(TOL)) {
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 1.0;
    int state = 1;
    BOOST_TEST(channel_model.pdf(observed, state) == 2.4912255069616864);
}

BOOST_AUTO_TEST_CASE(pdf_state_eq_obs_ne_one_test, *tolerance(TOL)) {
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 1.3;
    int state = 1;
    BOOST_TEST(channel_model.pdf(observed, state) == 0.43085303703574312);
}

BOOST_AUTO_TEST_CASE(sigma_test, *tolerance(TOL)) {
    ChannelModel channel_model;
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    BOOST_TEST(channel_model.sigma(3)
               == sqrt(0.00667 * 0.00667 + 3.0 * .16 * .16));
}

BOOST_AUTO_TEST_CASE(distance_p_bleach_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.bg_sig = 0.00667;
    cm1.mu = 0.5;
    cm1.sig = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.66;
    cm2.p_dud = 0.5;
    cm2.bg_sig = 0.00667;
    cm2.mu = 0.5;
    cm2.sig = 0.5;

    BOOST_TEST(cm1.distance(cm2) == (0.66 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_p_dud_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.bg_sig = 0.00667;
    cm1.mu = 0.5;
    cm1.sig = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.66;
    cm2.bg_sig = 0.00667;
    cm2.mu = 0.5;
    cm2.sig = 0.5;

    BOOST_TEST(cm1.distance(cm2) == (0.66 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_bg_sig_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.bg_sig = 0.5;
    cm1.mu = 0.5;
    cm1.sig = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.5;
    cm2.bg_sig = 0.66;
    cm2.mu = 0.5;
    cm2.sig = 0.5;

    BOOST_TEST(cm1.distance(cm2) == (0.66 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_mu_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.bg_sig = 0.00667;
    cm1.mu = 0.5;
    cm1.sig = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.5;
    cm2.bg_sig = 0.00667;
    cm2.mu = 0.66;
    cm2.sig = 0.5;

    BOOST_TEST(cm1.distance(cm2) == (0.66 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_sig_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.bg_sig = 0.00667;
    cm1.mu = 0.5;
    cm1.sig = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.5;
    cm2.bg_sig = 0.00667;
    cm2.mu = 0.5;
    cm2.sig = 0.66;

    BOOST_TEST(cm1.distance(cm2) == (0.66 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_max_no_sum_test, *tolerance(TOL)) {
    ChannelModel cm1;
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;
    cm1.bg_sig = 0.00667;
    cm1.mu = 0.66;
    cm1.sig = 0.5;

    ChannelModel cm2;
    cm2.p_bleach = 0.7;
    cm2.p_dud = 0.7;
    cm2.bg_sig = 0.00667;
    cm2.mu = 0.66;
    cm2.sig = 0.7;

    BOOST_TEST(cm1.distance(cm2) == (0.7 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.7 - 0.5));
}

BOOST_AUTO_TEST_SUITE_END()  // channel_model_suite
BOOST_AUTO_TEST_SUITE_END()  // model_suite
BOOST_AUTO_TEST_SUITE_END()  // parameterization_suite

}  // namespace whatprot
