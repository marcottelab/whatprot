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
#include "stuck-dye-state-vector.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(state_vector_suite)
BOOST_AUTO_TEST_SUITE(stuck_dye_state_vector_suite)

BOOST_AUTO_TEST_CASE(initialize_from_start_test) {
    StuckDyeStateVector sdsv;
    sdsv.initialize_from_start();
    BOOST_TEST(sdsv.dye == 1.0);
    BOOST_TEST(sdsv.no_dye == 0.0);
}

BOOST_AUTO_TEST_CASE(initialize_from_finish_test) {
    StuckDyeStateVector sdsv;
    sdsv.initialize_from_finish();
    BOOST_TEST(sdsv.dye == 1.0);
    BOOST_TEST(sdsv.no_dye == 1.0);
}

BOOST_AUTO_TEST_CASE(sum_test) {
    StuckDyeStateVector sdsv;
    sdsv.dye = 2.0;
    sdsv.no_dye = 2.0;
    BOOST_TEST(sdsv.sum() == 4.0);
}

BOOST_AUTO_TEST_CASE(source_test) {
    StuckDyeStateVector sdsv;
    sdsv.dye = 3.14;
    sdsv.no_dye = 2.12;
    BOOST_TEST(sdsv.source() == 3.14);
}

BOOST_AUTO_TEST_SUITE_END()  // stuck_dye_state_vector_suite
BOOST_AUTO_TEST_SUITE_END()  // state_vector_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
