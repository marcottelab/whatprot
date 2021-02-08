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
#include "start.h"

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(start_suite)

BOOST_AUTO_TEST_CASE(trivial_test, *tolerance(TOL)) {
    Start start;
    int order = 1;
    int* shape = new int[order];
    shape[0] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    tsr[loc] = -1000.0;
    int edmans;  // should be ignored.
    start.forward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 1.0);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(many_timesteps_test, *tolerance(TOL)) {
    Start start;
    int order = 1;
    int* shape = new int[order];
    shape[0] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    tsr[loc] = -1000.0;
    int edmans;  // should be ignored.
    start.forward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 1.0);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(many_dye_counts_test, *tolerance(TOL)) {
    Start start;
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = -1000.0;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = -1000.0;  // loc is {0, 2}
    int edmans;  // should be ignored.
    start.forward(tsr, &edmans, &tsr);
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 1.0);  // loc is {0, 2}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(more_dye_colors_test, *tolerance(TOL)) {
    Start start;
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {0, 1, 1}
    int edmans;  // should be ignored.
    start.forward(tsr, &edmans, &tsr);
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 1.0);  // loc is {0, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(everything_together_test, *tolerance(TOL)) {
    Start start;
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {0, 1, 1}
    int edmans;  // should be ignored.
    start.forward(tsr, &edmans, &tsr);
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 1.0);  // loc is {0, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_SUITE_END()  // start_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
