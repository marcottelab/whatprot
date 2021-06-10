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
#include "common/radiometry.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/step.h"
#include "hmm/step/stuck-dye-emission.h"
#include "hmm/step/stuck-dye-transition.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"

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
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.01;
    seq_model.p_detach = 0.02;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.03;
        seq_model.channel_models[i]->p_dud = 0.04;
        seq_model.channel_models[i]->mu = log(1.0);
        seq_model.channel_models[i]->sigma = 0.05;
        seq_model.channel_models[i]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[i]->p_stuck_dye_loss = 0.08;
    }
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(seq_model, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 4;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 1.0;
    r(0, 1) = 1.0;
    r(1, 0) = 1.0;
    r(1, 1) = 1.0;
    r(2, 0) = 1.0;
    r(2, 1) = 1.0;
    r(3, 0) = 1.0;
    r(3, 1) = 1.0;
    RadiometryPrecomputations radiometry_precomputations(
            r, seq_model, max_num_dyes);
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
    BOOST_TEST(*step
               == &universal_precomputations.stuck_dye_transitions[channel]);
    step++;
    BOOST_TEST(*step
               == &radiometry_precomputations.stuck_dye_emissions[channel]);
    step++;
    BOOST_TEST(*step
               == &universal_precomputations.stuck_dye_transitions[channel]);
    step++;
    BOOST_TEST(*step
               == &radiometry_precomputations.stuck_dye_emissions[channel]);
    step++;
    BOOST_TEST(*step
               == &universal_precomputations.stuck_dye_transitions[channel]);
    step++;
    BOOST_TEST(*step
               == &radiometry_precomputations.stuck_dye_emissions[channel]);
}

BOOST_AUTO_TEST_CASE(probability_more_involved_test, *tolerance(TOL)) {
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on April 7, 2021.
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.01;
    seq_model.p_detach = 0.02;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.03;
        seq_model.channel_models[i]->p_dud = 0.04;
        seq_model.channel_models[i]->mu = log(1.0);
        seq_model.channel_models[i]->sigma = 0.05;
        seq_model.channel_models[i]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[i]->p_stuck_dye_loss = 0.08;
    }
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(seq_model, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 3;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 1.0;
    r(0, 1) = 0.0;
    r(1, 0) = 1.0;
    r(1, 1) = 0.0;
    r(2, 0) = 0.0;
    r(2, 1) = 0.0;
    RadiometryPrecomputations radiometry_precomputations(
            r, seq_model, max_num_dyes);
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
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.01;
    seq_model.p_detach = 0.02;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.03;
        seq_model.channel_models[i]->p_dud = 0.04;
        seq_model.channel_models[i]->mu = log(1.0);
        seq_model.channel_models[i]->sigma = 0.05;
        seq_model.channel_models[i]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[i]->p_stuck_dye_loss = 0.08;
    }
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(seq_model, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 3;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 1.0;
    r(0, 1) = 0.0;
    r(1, 0) = 1.0;
    r(1, 1) = 0.0;
    r(2, 0) = 0.0;
    r(2, 1) = 0.0;
    RadiometryPrecomputations radiometry_precomputations(
            r, seq_model, max_num_dyes);
    int channel = 0;
    StuckDyeHMM hmm(num_timesteps,
                    num_channels,
                    channel,
                    radiometry_precomputations,
                    universal_precomputations);
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    smf.channel_fits.push_back(new ChannelModelFitter());
    hmm.improve_fit(&smf);
    // There are no BOOST_TEST statements because setting up a proper test for
    // this function is very difficult. We still have the test though as a no
    // crash test.
}

BOOST_AUTO_TEST_SUITE_END()  // stuck_dye_hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
