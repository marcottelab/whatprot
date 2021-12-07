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
#include "sequencing-model.h"

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
BOOST_AUTO_TEST_SUITE(sequencing_model_suite)

BOOST_AUTO_TEST_CASE(relative_distance_p_edman_failure_test, *tolerance(TOL)) {
    SequencingModel sm1;
    sm1.p_edman_failure = 0.5;
    sm1.p_detach = 0.5;

    SequencingModel sm2;
    sm2.p_edman_failure = 0.66;
    sm2.p_detach = 0.5;

    BOOST_TEST(sm1.relative_distance(sm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(sm2.relative_distance(sm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_p_detach_test, *tolerance(TOL)) {
    SequencingModel sm1;
    sm1.p_edman_failure = 0.5;
    sm1.p_detach = 0.5;

    SequencingModel sm2;
    sm2.p_edman_failure = 0.5;
    sm2.p_detach = 0.66;

    BOOST_TEST(sm1.relative_distance(sm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(sm2.relative_distance(sm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_with_channel_model_test,
                     *tolerance(TOL)) {
    SequencingModel sm1;
    sm1.p_edman_failure = 0.5;
    sm1.p_detach = 0.5;
    sm1.channel_models.push_back(new ChannelModel());
    sm1.channel_models[0]->p_bleach = 0.5;
    sm1.channel_models[0]->p_dud = 0.5;
    sm1.channel_models[0]->mu = 0.5;
    sm1.channel_models[0]->sig = 0.5;
    sm1.channel_models[0]->stuck_dye_ratio = 0.5;
    sm1.channel_models[0]->p_stuck_dye_loss = 0.5;

    SequencingModel sm2;
    sm2.p_edman_failure = 0.5;
    sm2.p_detach = 0.5;
    sm2.channel_models.push_back(new ChannelModel());
    sm2.channel_models[0]->p_bleach = 0.5;
    sm2.channel_models[0]->p_dud = 0.66;
    sm2.channel_models[0]->mu = 0.5;
    sm2.channel_models[0]->sig = 0.5;
    sm2.channel_models[0]->stuck_dye_ratio = 0.5;
    sm2.channel_models[0]->p_stuck_dye_loss = 0.5;

    BOOST_TEST(sm1.relative_distance(sm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(sm2.relative_distance(sm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_with_two_channel_models_test,
                     *tolerance(TOL)) {
    SequencingModel sm1;
    sm1.p_edman_failure = 0.5;
    sm1.p_detach = 0.5;
    sm1.channel_models.push_back(new ChannelModel());
    sm1.channel_models[0]->p_bleach = 0.5;
    sm1.channel_models[0]->p_dud = 0.5;
    sm1.channel_models[0]->mu = 0.5;
    sm1.channel_models[0]->sig = 0.5;
    sm1.channel_models[0]->stuck_dye_ratio = 0.5;
    sm1.channel_models[0]->p_stuck_dye_loss = 0.5;
    sm1.channel_models.push_back(new ChannelModel());
    sm1.channel_models[1]->p_bleach = 0.5;
    sm1.channel_models[1]->p_dud = 0.5;
    sm1.channel_models[1]->mu = 0.5;
    sm1.channel_models[1]->sig = 0.5;
    sm1.channel_models[1]->stuck_dye_ratio = 0.5;
    sm1.channel_models[1]->p_stuck_dye_loss = 0.5;

    SequencingModel sm2;
    sm2.p_edman_failure = 0.5;
    sm2.p_detach = 0.5;
    sm2.channel_models.push_back(new ChannelModel());
    sm2.channel_models[0]->p_bleach = 0.5;
    sm2.channel_models[0]->p_dud = 0.5;
    sm2.channel_models[0]->mu = 0.5;
    sm2.channel_models[0]->sig = 0.5;
    sm2.channel_models[0]->stuck_dye_ratio = 0.5;
    sm2.channel_models[0]->p_stuck_dye_loss = 0.5;
    sm2.channel_models.push_back(new ChannelModel());
    sm2.channel_models[1]->p_bleach = 0.5;
    sm2.channel_models[1]->p_dud = 0.66;
    sm2.channel_models[1]->mu = 0.5;
    sm2.channel_models[1]->sig = 0.5;
    sm2.channel_models[1]->stuck_dye_ratio = 0.5;
    sm2.channel_models[1]->p_stuck_dye_loss = 0.5;

    BOOST_TEST(sm1.relative_distance(sm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(sm2.relative_distance(sm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_CASE(relative_distance_max_no_sum_test, *tolerance(TOL)) {
    SequencingModel sm1;
    sm1.p_edman_failure = 0.5;
    sm1.p_detach = 0.5;
    sm1.channel_models.push_back(new ChannelModel());
    sm1.channel_models[0]->p_bleach = 0.5;
    sm1.channel_models[0]->p_dud = 0.5;
    sm1.channel_models[0]->mu = 0.5;
    sm1.channel_models[0]->sig = 0.5;
    sm1.channel_models[0]->stuck_dye_ratio = 0.5;
    sm1.channel_models[0]->p_stuck_dye_loss = 0.5;
    sm1.channel_models.push_back(new ChannelModel());
    sm1.channel_models[1]->p_bleach = 0.5;
    sm1.channel_models[1]->p_dud = 0.5;
    sm1.channel_models[1]->mu = 0.5;
    sm1.channel_models[1]->sig = 0.5;
    sm1.channel_models[1]->stuck_dye_ratio = 0.5;
    sm1.channel_models[1]->p_stuck_dye_loss = 0.5;

    SequencingModel sm2;
    sm2.p_edman_failure = 0.66;
    sm2.p_detach = 0.66;
    sm2.channel_models.push_back(new ChannelModel());
    sm2.channel_models[0]->p_bleach = 0.5;
    sm2.channel_models[0]->p_dud = 0.5;
    sm2.channel_models[0]->mu = 0.5;
    sm2.channel_models[0]->sig = 0.5;
    sm2.channel_models[0]->stuck_dye_ratio = 0.66;
    sm2.channel_models[0]->p_stuck_dye_loss = 0.5;
    sm2.channel_models.push_back(new ChannelModel());
    sm2.channel_models[1]->p_bleach = 0.5;
    sm2.channel_models[1]->p_dud = 0.66;
    sm2.channel_models[1]->mu = 0.5;
    sm2.channel_models[1]->sig = 0.5;
    sm2.channel_models[1]->stuck_dye_ratio = 0.5;
    sm2.channel_models[1]->p_stuck_dye_loss = 0.5;

    BOOST_TEST(sm1.relative_distance(sm2) == (0.66 - 0.5) / 0.5);
    BOOST_TEST(sm2.relative_distance(sm1) == (0.66 - 0.5) / 0.66);
}

BOOST_AUTO_TEST_SUITE_END()  // sequencing_model_suite
BOOST_AUTO_TEST_SUITE_END()  // model_suite
BOOST_AUTO_TEST_SUITE_END()  // parameterization_suite

}  // namespace whatprot
