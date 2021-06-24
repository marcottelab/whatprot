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
#include "detach-transition.h"
#include "hmm/state-vector/peptide-state-vector.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(detach_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    BOOST_CHECK_EQUAL(dt.p_detach, p_detach);
}

BOOST_AUTO_TEST_CASE(forward_in_place_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 1.0;
    unsigned int edmans = 0;
    dt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 1.0);
}

BOOST_AUTO_TEST_CASE(forward_in_place_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.3;
    psv.tensor[{0, 1}] = 0.7;
    unsigned int edmans = 0;
    dt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.3 + 0.7 * p_detach);
    BOOST_TEST((psv.tensor[{0, 1}]) == 0.7 * (1 - p_detach));
}

BOOST_AUTO_TEST_CASE(forward_in_place_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.3;
    psv.tensor[{0, 1}] = 0.6;
    psv.tensor[{0, 2}] = 0.1;
    unsigned int edmans = 0;
    dt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.3 + (0.6 + 0.1) * p_detach);
    BOOST_TEST((psv.tensor[{0, 1}]) == 0.6 * (1 - p_detach));
    BOOST_TEST((psv.tensor[{0, 2}]) == 0.1 * (1 - p_detach));
}

BOOST_AUTO_TEST_CASE(forward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.1;
    psv.tensor[{0, 1}] = 0.2;
    psv.tensor[{1, 0}] = 0.3;
    psv.tensor[{1, 1}] = 0.4;
    psv.tensor[{2, 0}] = 0.5;
    psv.tensor[{2, 1}] = 0.6;
    unsigned int edmans = 2;
    dt.forward(&edmans, &psv);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    BOOST_TEST((psv.tensor[{0, 1}]) == 0.2 * (1 - p_detach));
    BOOST_TEST((psv.tensor[{1, 1}]) == 0.4 * (1 - p_detach));
    BOOST_TEST((psv.tensor[{2, 1}]) == 0.6 * (1 - p_detach));
    // Distribution between empty states is of no importance. We only care about
    // the sum.
    double sum_empties = 0.0;
    sum_empties += psv.tensor[{0, 0}];
    sum_empties += psv.tensor[{1, 0}];
    sum_empties += psv.tensor[{2, 0}];
    BOOST_TEST(sum_empties == 0.1 + 0.3 + 0.5 + (0.2 + 0.4 + 0.6) * p_detach);
}

BOOST_AUTO_TEST_CASE(forward_in_place_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0, 0}] = 0.1;
    psv.tensor[{0, 0, 1}] = 0.2;
    psv.tensor[{0, 1, 0}] = 0.3;
    psv.tensor[{0, 1, 1}] = 0.4;
    unsigned int edmans = 0;
    dt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0, 0}]) == 0.1 + (0.2 + 0.3 + 0.4) * p_detach);
    BOOST_TEST((psv.tensor[{0, 0, 1}]) == 0.2 * (1 - p_detach));
    BOOST_TEST((psv.tensor[{0, 1, 0}]) == 0.3 * (1 - p_detach));
    BOOST_TEST((psv.tensor[{0, 1, 1}]) == 0.4 * (1 - p_detach));
}

BOOST_AUTO_TEST_CASE(backward_in_place_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 1.0;
    unsigned int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 1.0);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 1.0;
    unsigned int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0}]) == 1.0);
}

BOOST_AUTO_TEST_CASE(backward_in_place_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.3;
    psv.tensor[{0, 1}] = 0.7;
    unsigned int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.3);
    BOOST_TEST((psv.tensor[{0, 1}]) == p_detach * 0.3 + (1 - p_detach) * 0.7);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    unsigned int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0}]) == 0.3);
    BOOST_TEST((psv2.tensor[{0, 1}]) == p_detach * 0.3 + (1 - p_detach) * 0.7);
}

BOOST_AUTO_TEST_CASE(backward_in_place_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.3;
    psv.tensor[{0, 1}] = 0.6;
    psv.tensor[{0, 2}] = 0.1;
    unsigned int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.3);
    BOOST_TEST((psv.tensor[{0, 1}]) == p_detach * 0.3 + (1 - p_detach) * 0.6);
    BOOST_TEST((psv.tensor[{0, 2}]) == p_detach * 0.3 + (1 - p_detach) * 0.1);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    unsigned int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0}]) == 0.3);
    BOOST_TEST((psv2.tensor[{0, 1}]) == p_detach * 0.3 + (1 - p_detach) * 0.6);
    BOOST_TEST((psv2.tensor[{0, 2}]) == p_detach * 0.3 + (1 - p_detach) * 0.1);
}

BOOST_AUTO_TEST_CASE(backward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.88;
    psv.tensor[{0, 1}] = 0.2;
    psv.tensor[{1, 0}] = 0.88;
    psv.tensor[{1, 1}] = 0.4;
    psv.tensor[{2, 0}] = 0.88;
    psv.tensor[{2, 1}] = 0.6;
    unsigned int edmans = 2;
    dt.backward(psv, &edmans, &psv);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.88);
    BOOST_TEST((psv.tensor[{0, 1}]) == p_detach * 0.88 + (1 - p_detach) * 0.2);
    BOOST_TEST((psv.tensor[{1, 0}]) == 0.88);
    BOOST_TEST((psv.tensor[{1, 1}]) == p_detach * 0.88 + (1 - p_detach) * 0.4);
    BOOST_TEST((psv.tensor[{2, 0}]) == 0.88);
    BOOST_TEST((psv.tensor[{2, 1}]) == p_detach * 0.88 + (1 - p_detach) * 0.6);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.88;
    psv1.tensor[{0, 1}] = 0.2;
    psv1.tensor[{1, 0}] = 0.88;
    psv1.tensor[{1, 1}] = 0.4;
    psv1.tensor[{2, 0}] = 0.88;
    psv1.tensor[{2, 1}] = 0.6;
    unsigned int edmans = 2;
    dt.backward(psv1, &edmans, &psv2);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    BOOST_TEST((psv2.tensor[{0, 0}]) == 0.88);
    BOOST_TEST((psv2.tensor[{0, 1}]) == p_detach * 0.88 + (1 - p_detach) * 0.2);
    BOOST_TEST((psv2.tensor[{1, 0}]) == 0.88);
    BOOST_TEST((psv2.tensor[{1, 1}]) == p_detach * 0.88 + (1 - p_detach) * 0.4);
    BOOST_TEST((psv2.tensor[{2, 0}]) == 0.88);
    BOOST_TEST((psv2.tensor[{2, 1}]) == p_detach * 0.88 + (1 - p_detach) * 0.6);
}

BOOST_AUTO_TEST_CASE(backward_in_place_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0, 0}] = 0.1;
    psv.tensor[{0, 0, 1}] = 0.2;
    psv.tensor[{0, 1, 0}] = 0.3;
    psv.tensor[{0, 1, 1}] = 0.4;
    unsigned int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0, 0}]) == 0.1);
    BOOST_TEST((psv.tensor[{0, 0, 1}])
               == p_detach * 0.1 + (1 - p_detach) * 0.2);
    BOOST_TEST((psv.tensor[{0, 1, 0}])
               == p_detach * 0.1 + (1 - p_detach) * 0.3);
    BOOST_TEST((psv.tensor[{0, 1, 1}])
               == p_detach * 0.1 + (1 - p_detach) * 0.4);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 1}] = 0.2;
    psv1.tensor[{0, 1, 0}] = 0.3;
    psv1.tensor[{0, 1, 1}] = 0.4;
    unsigned int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0, 0}]) == 0.1);
    BOOST_TEST((psv2.tensor[{0, 0, 1}])
               == p_detach * 0.1 + (1 - p_detach) * 0.2);
    BOOST_TEST((psv2.tensor[{0, 1, 0}])
               == p_detach * 0.1 + (1 - p_detach) * 0.3);
    BOOST_TEST((psv2.tensor[{0, 1, 1}])
               == p_detach * 0.1 + (1 - p_detach) * 0.4);
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    fpsv.tensor[{0, 2}] = 0.91;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    bpsv.tensor[{0, 2}] = 0.92;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    nbpsv.tensor[{0, 2}] = 0.93;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 0.31 * 0.32 + 0.71 * 0.72 + 0.91 * 0.92;
    SequencingModelFitter smf;
    dt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    BOOST_TEST(smf.p_detach_fit.get()
               == (0.71 * p_detach * 0.33 + 0.91 * p_detach * 0.33)
                          / (0.71 * 0.72 + 0.91 * 0.92));
}

BOOST_AUTO_TEST_CASE(improve_fit_twice_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector fpsv1(order, shape);
    fpsv1.tensor[{0, 0}] = 0.31;
    fpsv1.tensor[{0, 1}] = 0.71;
    fpsv1.tensor[{0, 2}] = 0.91;
    PeptideStateVector bpsv1(order, shape);
    bpsv1.tensor[{0, 0}] = 0.32;
    bpsv1.tensor[{0, 1}] = 0.72;
    bpsv1.tensor[{0, 2}] = 0.92;
    PeptideStateVector nbpsv1(order, shape);
    nbpsv1.tensor[{0, 0}] = 0.33;
    nbpsv1.tensor[{0, 1}] = 0.73;
    nbpsv1.tensor[{0, 2}] = 0.93;
    PeptideStateVector fpsv2(order, shape);
    fpsv2.tensor[{0, 0}] = 0.231;
    fpsv2.tensor[{0, 1}] = 0.271;
    fpsv2.tensor[{0, 2}] = 0.291;
    PeptideStateVector bpsv2(order, shape);
    bpsv2.tensor[{0, 0}] = 0.232;
    bpsv2.tensor[{0, 1}] = 0.272;
    bpsv2.tensor[{0, 2}] = 0.292;
    PeptideStateVector nbpsv2(order, shape);
    nbpsv2.tensor[{0, 0}] = 0.233;
    nbpsv2.tensor[{0, 1}] = 0.273;
    nbpsv2.tensor[{0, 2}] = 0.293;
    delete[] shape;
    unsigned int edmans = 0;
    double prob1 = 0.31 * 0.32 + 0.71 * 0.72 + 0.91 * 0.92;
    double prob2 = 0.231 * 0.232 + 0.271 * 0.272 + 0.291 * 0.292;
    SequencingModelFitter smf;
    dt.improve_fit(fpsv1, bpsv1, nbpsv1, edmans, prob1, &smf);
    dt.improve_fit(fpsv2, bpsv2, nbpsv2, edmans, prob2, &smf);
    BOOST_TEST(
            smf.p_detach_fit.get()
            == ((0.71 * p_detach * 0.33 + 0.91 * p_detach * 0.33) / prob1
                + (0.271 * p_detach * 0.233 + 0.291 * p_detach * 0.233) / prob2)
                       / ((0.71 * 0.72 + 0.91 * 0.92) / prob1
                          + (0.271 * 0.272 + 0.291 * 0.292) / prob2));
}

BOOST_AUTO_TEST_SUITE_END()  // detach_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
