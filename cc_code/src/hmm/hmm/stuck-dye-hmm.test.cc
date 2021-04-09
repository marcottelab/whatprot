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
#include "stuck-dye-hmm.h"

// Standard C++ library headers:
#include <cmath>
#include <vector>

// Local project headers:
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/fit/error-model-fitter.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/step.h"
#include "hmm/step/stuck-dye-emission.h"
#include "hmm/step/stuck-dye-transition.h"

namespace {
using boost::unit_test::tolerance;
using std::log;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

namespace whatprot {

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(stuck_dye_hmm_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double p_edman_failure = 0.01;
    double p_detach = 0.02;
    double p_bleach = 0.03;
    double p_dud = 0.04;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = 0.05;
    double stuck_dye_ratio = 0.5;
    double p_stuck_dye_loss = 0.08;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio,
                  p_stuck_dye_loss);
    int num_channels = 2;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(em, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 4;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 1.0;
    r(0, 1) = 1.0;
    r(1, 0) = 1.0;
    r(1, 1) = 1.0;
    r(2, 0) = 1.0;
    r(2, 1) = 1.0;
    r(3, 0) = 1.0;
    r(3, 1) = 1.0;
    RadiometryPrecomputations radiometry_precomputations(r, em, max_num_dyes);
    int channel = 0;
    StuckDyeHMM hmm(num_timesteps,
                    num_channels,
                    channel,
                    radiometry_precomputations,
                    universal_precomputations);
    BOOST_ASSERT(hmm.steps.size() == 2 * num_timesteps - 1);
    vector<const Step<StuckDyeStateVector>*>::iterator step = hmm.steps.begin();
    BOOST_TEST(*step
               == &radiometry_precomputations.stuck_dye_emissions[channel]);
    step++;
    BOOST_TEST(*step == &universal_precomputations.stuck_dye_transition);
    step++;
    BOOST_TEST(*step
               == &radiometry_precomputations.stuck_dye_emissions[channel]);
    step++;
    BOOST_TEST(*step == &universal_precomputations.stuck_dye_transition);
    step++;
    BOOST_TEST(*step
               == &radiometry_precomputations.stuck_dye_emissions[channel]);
    step++;
    BOOST_TEST(*step == &universal_precomputations.stuck_dye_transition);
    step++;
    BOOST_TEST(*step
               == &radiometry_precomputations.stuck_dye_emissions[channel]);
}

BOOST_AUTO_TEST_CASE(probability_trivial_test, *tolerance(TOL)) {
    double p_edman_failure = 0.01;
    double p_detach = 0.02;
    double p_bleach = 0.03;
    double p_dud = 0.04;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = 0.05;
    double stuck_dye_ratio = 0.5;
    double p_stuck_dye_loss = 0.08;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio,
                  p_stuck_dye_loss);
    int num_channels = 2;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(em, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 1;
    Radiometry r(num_timesteps, num_channels);
    RadiometryPrecomputations radiometry_precomputations(r, em, max_num_dyes);
    int channel = 0;
    StuckDyeHMM hmm(num_timesteps,
                    num_channels,
                    channel,
                    radiometry_precomputations,
                    universal_precomputations);
    BOOST_TEST(hmm.probability() == 1.0);
}

BOOST_AUTO_TEST_CASE(probability_sum_to_one_test, *tolerance(TOL)) {
    // The idea here is that if the emission just causes a multiplication by one
    // no matter what, then the transitions shouldn't matter; they should just
    // be reallocating a probability of one between different states.
    double p_edman_failure = 0.01;
    double p_detach = 0.02;
    double p_bleach = 0.03;
    double p_dud = 0.04;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = 0.05;
    double stuck_dye_ratio = 0.5;
    double p_stuck_dye_loss = 0.08;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio,
                  p_stuck_dye_loss);
    int num_channels = 2;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(em, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 3;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 0.0;
    r(0, 1) = 0.1;
    r(1, 0) = 1.0;
    r(1, 1) = 1.1;
    r(2, 0) = 2.0;
    r(2, 1) = 2.1;
    RadiometryPrecomputations radiometry_precomputations(r, em, max_num_dyes);
    int channel = 0;
    StuckDyeHMM hmm(num_timesteps,
                    num_channels,
                    channel,
                    radiometry_precomputations,
                    universal_precomputations);
    BOOST_TEST(hmm.probability() == 1.0);
}

BOOST_AUTO_TEST_CASE(probability_more_involved_test, *tolerance(TOL)) {
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on April 7, 2021.
    double p_edman_failure = 0.01;
    double p_detach = 0.02;
    double p_bleach = 0.03;
    double p_dud = 0.04;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = 0.05;
    double stuck_dye_ratio = 0.5;
    double p_stuck_dye_loss = 0.08;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio,
                  p_stuck_dye_loss);
    int num_channels = 2;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(em, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 3;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 1.0;
    r(0, 1) = 0.0;
    r(1, 0) = 1.0;
    r(1, 1) = 0.0;
    r(2, 0) = 0.0;
    r(2, 1) = 0.0;
    RadiometryPrecomputations radiometry_precomputations(r, em, max_num_dyes);
    int channel = 0;
    StuckDyeHMM hmm(num_timesteps,
                    num_channels,
                    channel,
                    radiometry_precomputations,
                    universal_precomputations);
    BOOST_TEST(hmm.probability() == 4.6855215246253996);
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on April 7, 2021.
    double p_edman_failure = 0.01;
    double p_detach = 0.02;
    double p_bleach = 0.03;
    double p_dud = 0.04;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = log(1.0);
    double sigma = 0.05;
    double stuck_dye_ratio = 0.5;
    double p_stuck_dye_loss = 0.08;
    ErrorModel em(p_edman_failure,
                  p_detach,
                  p_bleach,
                  p_dud,
                  dist_type,
                  mu,
                  sigma,
                  stuck_dye_ratio,
                  p_stuck_dye_loss);
    int num_channels = 2;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(em, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 3;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 1.0;
    r(0, 1) = 0.0;
    r(1, 0) = 1.0;
    r(1, 1) = 0.0;
    r(2, 0) = 0.0;
    r(2, 1) = 0.0;
    RadiometryPrecomputations radiometry_precomputations(r, em, max_num_dyes);
    int channel = 0;
    StuckDyeHMM hmm(num_timesteps,
                    num_channels,
                    channel,
                    radiometry_precomputations,
                    universal_precomputations);
    ErrorModelFitter emf;
    hmm.improve_fit(&emf);
    // There are no BOOST_TEST statements because setting up a proper test for
    // this function is very difficult. We still have the test though as a no
    // crash test.
}

BOOST_AUTO_TEST_SUITE_END()  // stuck_dye_hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
