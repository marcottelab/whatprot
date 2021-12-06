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
#include "peptide-hmm.h"

// Standard C++ library headers:
#include <cmath>
#include <functional>
#include <typeinfo>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/binomial-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/step/edman-transition.h"
#include "hmm/step/peptide-emission.h"
#include "hmm/step/peptide-step.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "tensor/tensor.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using std::exp;
using std::function;
using std::log;
using std::sqrt;
using std::vector;
const double PI = 3.141592653589793238;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(peptide_hmm_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.01;
    seq_model.p_detach = 0.02;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.03;
        seq_model.channel_models[i]->p_dud = 0.04;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sigma = 0.05;
        seq_model.channel_models[i]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[i]->p_stuck_dye_loss = 0.08;
    }
    SequencingSettings seq_settings;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(seq_model, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 4;
    DyeSeq ds(num_channels, ".1.0.1.0.1");  // two in ch 0, three in ch 1.
    DyeSeqPrecomputations dye_seq_precomputations(
            ds, seq_model, num_timesteps, num_channels);
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
            r, seq_model, seq_settings, max_num_dyes);
    PeptideHMM hmm(num_timesteps,
                   num_channels,
                   dye_seq_precomputations,
                   radiometry_precomputations,
                   universal_precomputations);
    BOOST_ASSERT(hmm.tensor_shape.size() == 1 + num_channels);
    BOOST_TEST(hmm.tensor_shape[0] == num_timesteps);
    BOOST_TEST(hmm.tensor_shape[1] == 2u + 1u);  // extra for 0 & num dyes.
    BOOST_TEST(hmm.tensor_shape[2] == 3u + 1u);  // extra for 0 & num dyes.
    BOOST_ASSERT(hmm.steps.size()
                 == (3 + num_channels) * (num_timesteps - 1) + 1
                            + num_channels);
    vector<PeptideStep*>::iterator step = hmm.steps.begin();
    BOOST_TEST(typeid(**step).name() == typeid(DudTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(DudTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(PeptideEmission).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(DetachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(BleachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(BleachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(EdmanTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(PeptideEmission).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(DetachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(BleachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(BleachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(EdmanTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(PeptideEmission).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(DetachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(BleachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(BleachTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(EdmanTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(PeptideEmission).name());
}

BOOST_AUTO_TEST_CASE(probability_more_involved_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.05;
        seq_model.channel_models[i]->p_dud = 0.07;
        seq_model.channel_models[i]->bg_sigma = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sigma = 0.16;
        seq_model.channel_models[i]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[i]->p_stuck_dye_loss = 0.08;
    }
    SequencingSettings seq_settings;
    int max_num_dyes = 5;
    UniversalPrecomputations up(seq_model, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 3;
    DyeSeq ds(num_channels, "10.01111");  // two in ch 0, five in ch 1.
    DyeSeqPrecomputations dsp(ds, seq_model, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 5.0;
    r(0, 1) = 2.0;
    r(1, 0) = 5.0;
    r(1, 1) = 1.0;
    r(2, 0) = 4.0;
    r(2, 1) = 1.0;
    RadiometryPrecomputations rp(r, seq_model, seq_settings, max_num_dyes);
    PeptideHMM hmm(num_timesteps, num_channels, dsp, rp, up);
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on September 8, 2021.
    BOOST_TEST(hmm.probability() == 1.876822091893613e-96);
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.01;
    seq_model.p_detach = 0.02;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.03;
        seq_model.channel_models[i]->p_dud = 0.04;
        seq_model.channel_models[i]->bg_sigma = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sigma = 0.05;
        seq_model.channel_models[i]->stuck_dye_ratio = 0.5;
        seq_model.channel_models[i]->p_stuck_dye_loss = 0.08;
    }
    SequencingSettings seq_settings;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(seq_model, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 4;
    DyeSeq ds(num_channels, ".1.0.1.0.1");  // two in ch 0, three in ch 1.
    DyeSeqPrecomputations dye_seq_precomputations(
            ds, seq_model, num_timesteps, num_channels);
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
            r, seq_model, seq_settings, max_num_dyes);
    PeptideHMM hmm(num_timesteps,
                   num_channels,
                   dye_seq_precomputations,
                   radiometry_precomputations,
                   universal_precomputations);
    SequencingModelFitter smf(num_channels);
    hmm.improve_fit(&smf);
    // There are no BOOST_TEST statements because setting up a proper test for
    // this function is very difficult. We still have the test though as a no
    // crash test.
}

BOOST_AUTO_TEST_SUITE_END()  // peptide_hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
