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
#include "dud-transition.h"

// Local project headers:
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(dud_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double q = 0.2;
    int channel = -1;  // can be ignored for this test.
    DudTransition dt(q, channel);
    BOOST_TEST(dt.q == q);
    BOOST_TEST(dt.channel == channel);
}

BOOST_AUTO_TEST_CASE(improve_fit_basic_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    DudTransition dt(q, channel);
    dt.reserve(1);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector fpsv(order, shape);
    PeptideStateVector bpsv(order, shape);
    PeptideStateVector nbpsv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    fpsv.tensor[loc] = 0.31;  // loc is {0, 0}
    bpsv.tensor[loc] = 0.32;
    nbpsv.tensor[loc] = 0.33;
    loc[1] = 1;
    fpsv.tensor[loc] = 0.71;  // loc is {0, 1}
    bpsv.tensor[loc] = 0.72;
    nbpsv.tensor[loc] = 0.73;
    delete[] loc;
    int edmans = 0;
    double probability = 1.0;
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    dt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    BOOST_TEST(smf.channel_fits[0]->p_dud_fit.get()
               == (0.71 * q * 0.33) / (0.71 * 0.72));
}

BOOST_AUTO_TEST_SUITE_END()  // dud_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
