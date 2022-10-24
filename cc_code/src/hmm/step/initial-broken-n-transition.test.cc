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
#include "initial-broken-n-transition.h"

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
BOOST_AUTO_TEST_SUITE(initial_broken_n_transition_suite)

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    double p_break_n = 0.07;
    InitialBrokenNTransition ibnt(p_break_n);
    ibnt.pruned_range.min = {0, 0};
    ibnt.pruned_range.max = {1, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    fpsv.broken_n_tensor[{0, 0}] = 0.031;
    fpsv.broken_n_tensor[{0, 1}] = 0.071;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    bpsv.broken_n_tensor[{0, 0}] = 0.032;
    bpsv.broken_n_tensor[{0, 1}] = 0.072;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    nbpsv.broken_n_tensor[{0, 0}] = 0.033;
    nbpsv.broken_n_tensor[{0, 1}] = 0.073;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 1.0;
    SequencingModelFitter smf;
    ibnt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    BOOST_TEST(smf.p_initial_break_n_fit.get()
               == (0.31 * p_break_n * 0.033 + 0.71 * p_break_n * 0.073)
                          / (0.31 * 0.32 + 0.71 * 0.72));
}

BOOST_AUTO_TEST_SUITE_END()  // initial_broken_n_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
