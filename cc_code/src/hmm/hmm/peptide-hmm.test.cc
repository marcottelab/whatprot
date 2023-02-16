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
#include <limits>
#include <typeinfo>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/bleach-transition.h"
#include "hmm/step/cyclic-broken-n-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/step/dud-transition.h"
#include "hmm/step/edman-transition.h"
#include "hmm/step/initial-broken-n-transition.h"
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
    seq_model.p_initial_break_n = 0.07;
    seq_model.p_cyclic_break_n = 0.025;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.03;
        seq_model.channel_models[i]->p_dud = 0.04;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sig = 0.05;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
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
    BOOST_ASSERT(hmm.steps.size()
                 == (4 + num_channels) * (num_timesteps - 1) + 2
                            + num_channels);
    vector<PeptideStep*>::iterator step = hmm.steps.begin();
    BOOST_TEST(typeid(**step).name()
               == typeid(InitialBrokenNTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(DudTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(DudTransition).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(PeptideEmission).name());
    step++;
    BOOST_TEST(typeid(**step).name() == typeid(CyclicBrokenNTransition).name());
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
    BOOST_TEST(typeid(**step).name() == typeid(CyclicBrokenNTransition).name());
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
    BOOST_TEST(typeid(**step).name() == typeid(CyclicBrokenNTransition).name());
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

BOOST_AUTO_TEST_CASE(probability_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    seq_model.p_initial_break_n = 0.07;
    seq_model.p_cyclic_break_n = 0.025;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.05;
        seq_model.channel_models[i]->p_dud = 0.07;
        seq_model.channel_models[i]->bg_sig = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sig = 0.16;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    int max_num_dyes = 5;
    UniversalPrecomputations up(seq_model, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 3;
    DyeSeq ds(num_channels, "10.01111");  // two in ch 0, five in ch 1.
    DyeSeqPrecomputations dsp(ds, seq_model, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 2.0;
    r(0, 1) = 5.0;
    r(1, 0) = 1.0;
    r(1, 1) = 5.0;
    r(2, 0) = 1.0;
    r(2, 1) = 4.0;
    RadiometryPrecomputations rp(r, seq_model, seq_settings, max_num_dyes);
    PeptideHMM hmm(num_timesteps, num_channels, dsp, rp, up);
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on October 25, 2022.
    BOOST_TEST(hmm.probability() == 0.042073244987065869);
}

BOOST_AUTO_TEST_CASE(probability_distribution_tails_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    seq_model.p_initial_break_n = 0.07;
    seq_model.p_cyclic_break_n = 0.025;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.05;
        seq_model.channel_models[i]->p_dud = 0.07;
        seq_model.channel_models[i]->bg_sig = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sig = 0.16;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
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
    // giving the correct result on October 25, 2022.
    BOOST_TEST(hmm.probability() == 2.832048536598378e-96);
}

BOOST_AUTO_TEST_CASE(probability_detachment_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    seq_model.p_initial_break_n = 0.07;
    seq_model.p_cyclic_break_n = 0.025;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.05;
        seq_model.channel_models[i]->p_dud = 0.07;
        seq_model.channel_models[i]->bg_sig = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sig = 0.16;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    int max_num_dyes = 5;
    UniversalPrecomputations up(seq_model, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 3;
    DyeSeq ds(num_channels, "10.01111");  // two in ch 0, five in ch 1.
    DyeSeqPrecomputations dsp(ds, seq_model, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 2.0;
    r(0, 1) = 5.0;
    r(1, 0) = 0.0;
    r(1, 1) = 0.0;
    r(2, 0) = 0.0;
    r(2, 1) = 0.0;
    RadiometryPrecomputations rp(r, seq_model, seq_settings, max_num_dyes);
    PeptideHMM hmm(num_timesteps, num_channels, dsp, rp, up);
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on February 16 2023.
    BOOST_TEST(hmm.probability() == 758907.43743397435);
}

BOOST_AUTO_TEST_CASE(probability_with_cutoff_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    seq_model.p_initial_break_n = 0.07;
    seq_model.p_cyclic_break_n = 0.025;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.05;
        seq_model.channel_models[i]->p_dud = 0.07;
        seq_model.channel_models[i]->bg_sig = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sig = 0.16;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = 5.0;
    int max_num_dyes = 5;
    UniversalPrecomputations up(seq_model, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 3;
    DyeSeq ds(num_channels, "10.01111");  // two in ch 0, five in ch 1.
    DyeSeqPrecomputations dsp(ds, seq_model, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 2.0;
    r(0, 1) = 5.0;
    r(1, 0) = 1.0;
    r(1, 1) = 5.0;
    r(2, 0) = 1.0;
    r(2, 1) = 4.0;
    RadiometryPrecomputations rp(r, seq_model, seq_settings, max_num_dyes);
    PeptideHMM hmm(num_timesteps, num_channels, dsp, rp, up);
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on October 25, 2022.
    BOOST_TEST(hmm.probability() == 0.042073244682031316);
}

BOOST_AUTO_TEST_CASE(probability_with_cutoff_zero_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.06;
    seq_model.p_detach = 0.05;
    seq_model.p_initial_break_n = 0.07;
    seq_model.p_cyclic_break_n = 0.025;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.05;
        seq_model.channel_models[i]->p_dud = 0.07;
        seq_model.channel_models[i]->bg_sig = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sig = 0.16;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = 5.0;
    int max_num_dyes = 5;
    UniversalPrecomputations up(seq_model, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    unsigned int num_timesteps = 3;
    DyeSeq ds(num_channels, "10.01111");  // two in ch 0, five in ch 1.
    DyeSeqPrecomputations dsp(ds, seq_model, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 2.0;
    r(0, 1) = 1.0;
    r(1, 0) = 1.0;
    r(1, 1) = 5.0;
    r(2, 0) = 1.0;
    r(2, 1) = 4.0;
    RadiometryPrecomputations rp(r, seq_model, seq_settings, max_num_dyes);
    PeptideHMM hmm(num_timesteps, num_channels, dsp, rp, up);
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on October 25, 2022.
    BOOST_TEST(hmm.probability() == 0.0);
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    unsigned int num_channels = 2;
    SequencingModel seq_model;
    seq_model.p_edman_failure = 0.01;
    seq_model.p_detach = 0.02;
    seq_model.p_initial_break_n = 0.07;
    seq_model.p_cyclic_break_n = 0.025;
    for (unsigned int i = 0; i < num_channels; i++) {
        seq_model.channel_models.push_back(new ChannelModel());
        seq_model.channel_models[i]->p_bleach = 0.03;
        seq_model.channel_models[i]->p_dud = 0.04;
        seq_model.channel_models[i]->bg_sig = 0.00667;
        seq_model.channel_models[i]->mu = 1.0;
        seq_model.channel_models[i]->sig = 0.05;
    }
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
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
