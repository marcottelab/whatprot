/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// File under test:
#include "decaying-rate-model.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(parameterization_suite)
BOOST_AUTO_TEST_SUITE(model_suite)
BOOST_AUTO_TEST_SUITE(decaying_rate_model_suite)

BOOST_AUTO_TEST_CASE(distance_base_test, *tolerance(TOL)) {
    DecayingRateModel drm1;
    drm1.base = 0.5;
    drm1.initial = 0.5;
    drm1.initial_decay = 0.5;

    DecayingRateModel drm2;
    drm2.base = 0.66;
    drm2.initial = 0.5;
    drm2.initial_decay = 0.5;

    BOOST_TEST(drm1.distance(drm2) == (0.66 - 0.5));
    BOOST_TEST(drm2.distance(drm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_initial_test, *tolerance(TOL)) {
    DecayingRateModel drm1;
    drm1.base = 0.5;
    drm1.initial = 0.5;
    drm1.initial_decay = 0.5;

    DecayingRateModel drm2;
    drm2.base = 0.5;
    drm2.initial = 0.66;
    drm2.initial_decay = 0.5;

    BOOST_TEST(drm1.distance(drm2) == (0.66 - 0.5));
    BOOST_TEST(drm2.distance(drm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_CASE(distance_initial_decay_test, *tolerance(TOL)) {
    DecayingRateModel drm1;
    drm1.base = 0.5;
    drm1.initial = 0.5;
    drm1.initial_decay = 0.5;

    DecayingRateModel drm2;
    drm2.base = 0.5;
    drm2.initial = 0.5;
    drm2.initial_decay = 0.66;

    BOOST_TEST(drm1.distance(drm2) == (0.66 - 0.5));
    BOOST_TEST(drm2.distance(drm1) == (0.66 - 0.5));
}

BOOST_AUTO_TEST_SUITE_END()  // channel_model_suite
BOOST_AUTO_TEST_SUITE_END()  // model_suite
BOOST_AUTO_TEST_SUITE_END()  // parameterization_suite

}  // namespace whatprot
