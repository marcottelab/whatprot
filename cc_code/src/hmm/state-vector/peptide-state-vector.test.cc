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
#include "peptide-state-vector.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(state_vector_suite)
BOOST_AUTO_TEST_SUITE(peptide_state_vector_suite)

BOOST_AUTO_TEST_CASE(constructor_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 5;
    shape[2] = 7;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    BOOST_TEST(psv.tensor.order == 3);
    BOOST_TEST(psv.tensor.shape[0] == 3);
    BOOST_TEST(psv.tensor.shape[1] == 5);
    BOOST_TEST(psv.tensor.shape[2] == 7);
}

BOOST_AUTO_TEST_CASE(initialize_from_start_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    shape[2] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.initialize_from_start();
    BOOST_TEST(psv.tensor.values[0] == 0.0);
    BOOST_TEST(psv.tensor.values[1] == 0.0);
    BOOST_TEST(psv.tensor.values[2] == 0.0);
    BOOST_TEST(psv.tensor.values[3] == 0.0);
    BOOST_TEST(psv.tensor.values[4] == 0.0);
    BOOST_TEST(psv.tensor.values[5] == 1.0);
    BOOST_TEST(psv.tensor.values[6] == 0.0);
    BOOST_TEST(psv.tensor.values[7] == 0.0);
    BOOST_TEST(psv.tensor.values[8] == 0.0);
    BOOST_TEST(psv.tensor.values[9] == 0.0);
    BOOST_TEST(psv.tensor.values[10] == 0.0);
    BOOST_TEST(psv.tensor.values[11] == 0.0);
}

BOOST_AUTO_TEST_CASE(initialize_from_finish_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    shape[2] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.initialize_from_finish();
    BOOST_TEST(psv.tensor.values[0] == 1.0);
    BOOST_TEST(psv.tensor.values[1] == 1.0);
    BOOST_TEST(psv.tensor.values[2] == 1.0);
    BOOST_TEST(psv.tensor.values[3] == 1.0);
    BOOST_TEST(psv.tensor.values[4] == 1.0);
    BOOST_TEST(psv.tensor.values[5] == 1.0);
    BOOST_TEST(psv.tensor.values[6] == 1.0);
    BOOST_TEST(psv.tensor.values[7] == 1.0);
    BOOST_TEST(psv.tensor.values[8] == 1.0);
    BOOST_TEST(psv.tensor.values[9] == 1.0);
    BOOST_TEST(psv.tensor.values[10] == 1.0);
    BOOST_TEST(psv.tensor.values[11] == 1.0);
}

BOOST_AUTO_TEST_CASE(source_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    shape[2] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor.values[5] = 1.23456789;
    BOOST_TEST(psv.source() == 1.23456789);
}

BOOST_AUTO_TEST_SUITE_END()  // peptide_state_vector_suite
BOOST_AUTO_TEST_SUITE_END()  // state_vector_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
