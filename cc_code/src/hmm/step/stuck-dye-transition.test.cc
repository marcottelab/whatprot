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
#include "stuck-dye-transition.h"

// Local project headers:
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(stuck_dye_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    int channel = 0;
    double loss_rate = 0.035;
    StuckDyeTransition sdt(loss_rate, channel);
    BOOST_TEST(sdt.loss_rate == loss_rate);
}

BOOST_AUTO_TEST_CASE(forward_test, *tolerance(TOL)) {
    int channel = 0;
    double loss_rate = 0.035;
    StuckDyeTransition sdt(loss_rate, channel);
    unsigned int num_edmans = 0;
    StuckDyeStateVector sdsv1;
    sdsv1.dye = 0.3;
    sdsv1.no_dye = 0.7;
    StuckDyeStateVector* sdsv2 = sdt.forward(sdsv1, &num_edmans);
    BOOST_TEST(sdsv2->dye == 0.3 * (1 - loss_rate));
    BOOST_TEST(sdsv2->no_dye == 0.7 + 0.3 * loss_rate);
    BOOST_TEST(num_edmans == 1u);
}

BOOST_AUTO_TEST_CASE(backward_test, *tolerance(TOL)) {
    int channel = 0;
    double loss_rate = 0.035;
    StuckDyeTransition sdt(loss_rate, channel);
    unsigned int num_edmans = 1;
    StuckDyeStateVector input;
    input.dye = 0.3;
    input.no_dye = 0.7;
    StuckDyeStateVector* output = sdt.backward(input, &num_edmans);
    BOOST_TEST(output->dye == loss_rate * 0.7 + (1 - loss_rate) * 0.3);
    BOOST_TEST(output->no_dye == 0.7);
    BOOST_TEST(num_edmans == 0u);
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    int channel = 0;
    double loss_rate = 0.035;
    StuckDyeTransition sdt(loss_rate, channel);
    unsigned int num_edmans = 1;
    double probability = 0.98765;
    StuckDyeStateVector forward_sdsv;
    forward_sdsv.dye = 0.3;
    forward_sdsv.no_dye = 0.7;
    StuckDyeStateVector backward_sdsv;
    backward_sdsv.dye = 0.4;
    backward_sdsv.no_dye = 0.6;
    StuckDyeStateVector next_backward_sdsv;
    next_backward_sdsv.dye = 0.2;
    next_backward_sdsv.no_dye = 0.8;
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    sdt.improve_fit(forward_sdsv,
                    backward_sdsv,
                    next_backward_sdsv,
                    num_edmans,
                    probability,
                    &smf);
    BOOST_TEST(smf.channel_fits[0]->p_stuck_dye_loss_fit.get()
               == (0.3 * loss_rate * 0.8) / (0.3 * 0.4));
}

BOOST_AUTO_TEST_CASE(improve_fit_different_channel_test, *tolerance(TOL)) {
    int channel = 1;
    double loss_rate = 0.035;
    StuckDyeTransition sdt(loss_rate, channel);
    unsigned int num_edmans = 1;
    double probability = 0.98765;
    StuckDyeStateVector forward_sdsv;
    forward_sdsv.dye = 0.3;
    forward_sdsv.no_dye = 0.7;
    StuckDyeStateVector backward_sdsv;
    backward_sdsv.dye = 0.4;
    backward_sdsv.no_dye = 0.6;
    StuckDyeStateVector next_backward_sdsv;
    next_backward_sdsv.dye = 0.2;
    next_backward_sdsv.no_dye = 0.8;
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    smf.channel_fits.push_back(new ChannelModelFitter());
    sdt.improve_fit(forward_sdsv,
                    backward_sdsv,
                    next_backward_sdsv,
                    num_edmans,
                    probability,
                    &smf);
    BOOST_TEST(smf.channel_fits[0]->p_stuck_dye_loss_fit.numerator == 0);
    BOOST_TEST(smf.channel_fits[0]->p_stuck_dye_loss_fit.denominator == 0);
    BOOST_TEST(smf.channel_fits[1]->p_stuck_dye_loss_fit.get()
               == (0.3 * loss_rate * 0.8) / (0.3 * 0.4));
}

BOOST_AUTO_TEST_SUITE_END()  // stuck_dye_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
