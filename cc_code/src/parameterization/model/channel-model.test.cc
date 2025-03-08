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
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 0.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 0;
    BOOST_TEST(channel_model.pdf(observed, &counts[0]) == 59.811436342043883);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_state_zero_obs_one_test) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 1.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 0;
    // This is a regression test from January 22, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0]) == 0);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_state_one_obs_zero_test) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 0.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 1;
    // This is a regression test from January 22, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0])
               == 8.488175272749065e-09);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_state_one_obs_one_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 1.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 1;
    // This is a regression test from January 22, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0]) == 2.4912255069616864);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_state_eq_obs_ne_one_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    double observed = 1.3;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 1;
    // This is a regression test from January 22, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0]) == 0.43085303703574312);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_quench_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    channel_model.interactions[0] = 0.95;
    channel_model.flat_interactions[0] = 0.96;
    double observed = 2.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 2;
    // This is a regression test from January 24, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0]) == 1.3248217904301611);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_multiple_quench_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    channel_model.interactions[0] = 0.95;
    channel_model.flat_interactions[0] = 0.96;
    double observed = 3.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 3;
    // This is a regression test from January 24, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0]) == 0.46273832091183492);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_fret_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    channel_model.interactions[1] = 0.97;
    channel_model.flat_interactions[1] = 0.98;
    double observed = 3.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 1;
    counts[1] = 1;
    // This is a regression test from January 24, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0])
               == 9.961776308698737e-38);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_multiple_fret_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    channel_model.interactions[1] = 0.97;
    channel_model.flat_interactions[1] = 0.98;
    double observed = 3.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 1;
    counts[1] = 2;
    // This is a regression test from January 24, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0])
               == 5.8804666576569527e-40);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_quench_and_fret_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    channel_model.interactions[0] = 0.95;
    channel_model.flat_interactions[0] = 0.96;
    channel_model.interactions[1] = 0.97;
    channel_model.flat_interactions[1] = 0.98;
    double observed = 3.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 2;
    counts[1] = 1;
    // This is a regression test from January 24, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0])
               == 2.7722216687510138e-08);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(pdf_multiple_quench_multiple_fret_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    channel_model.interactions[0] = 0.95;
    channel_model.flat_interactions[0] = 0.96;
    channel_model.interactions[1] = 0.97;
    channel_model.flat_interactions[1] = 0.98;
    double observed = 3.0;
    unsigned int* counts = new unsigned int[num_channels];
    counts[0] = 3;
    counts[1] = 2;
    // This is a regression test from January 24, 2024.
    BOOST_TEST(channel_model.pdf(observed, &counts[0]) == 0.08307899998347866);
    delete[] counts;
}

BOOST_AUTO_TEST_CASE(sigma_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;
    ChannelModel channel_model(channel, num_channels);
    channel_model.p_bleach = .05;
    channel_model.p_dud = .10;
    channel_model.bg_sig = 0.00667;
    channel_model.mu = 1.0;
    channel_model.sig = .16;
    BOOST_TEST(channel_model.sigma(3)
               == sqrt(0.00667 * 0.00667 + 3.0 * .16 * .16));
}

BOOST_AUTO_TEST_CASE(distance_p_bleach_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;

    ChannelModel cm1(channel, num_channels);
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;

    ChannelModel cm2(channel, num_channels);
    cm2.p_bleach = 0.66;
    cm2.p_dud = 0.5;

    BOOST_TEST(cm1.distance(cm2) == (0.66 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_p_dud_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;

    ChannelModel cm1(channel, num_channels);
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;

    ChannelModel cm2(channel, num_channels);
    cm2.p_bleach = 0.5;
    cm2.p_dud = 0.66;

    BOOST_TEST(cm1.distance(cm2) == (0.66 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_max_no_sum_test, *tolerance(TOL)) {
    unsigned int num_channels = 1;
    unsigned int channel = 0;

    ChannelModel cm1(channel, num_channels);
    cm1.p_bleach = 0.5;
    cm1.p_dud = 0.5;

    ChannelModel cm2(channel, num_channels);
    cm2.p_bleach = 0.7;
    cm2.p_dud = 0.7;

    BOOST_TEST(cm1.distance(cm2) == (0.7 - 0.5));
    BOOST_TEST(cm2.distance(cm1) == (0.7 - 0.5));
}

BOOST_AUTO_TEST_SUITE_END()  // channel_model_suite
BOOST_AUTO_TEST_SUITE_END()  // model_suite
BOOST_AUTO_TEST_SUITE_END()  // parameterization_suite

}  // namespace whatprot
