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
#include "summation.h"

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(fwd_alg_suite)
BOOST_AUTO_TEST_SUITE(summation_suite)

BOOST_AUTO_TEST_CASE(trivial_test, *tolerance(TOL)) {
    Summation sum;
    int order = 1;
    int* shape = new int[order];
    shape[0] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    tsr[loc] = 3.14;
    int timestep = 0;
    BOOST_TEST(sum(&tsr, timestep) == 3.14);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(bigger_size_test, *tolerance(TOL)) {
    Summation sum;
    int order = 1;
    int* shape = new int[order];
    shape[0] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    tsr[loc] = 7.0;
    loc[0] = 1;
    tsr[loc] = 7.1;
    loc[0] = 2;
    tsr[loc] = 7.2;
    int timestep = 2;
    BOOST_TEST(sum(&tsr, timestep) == 7.0 + 7.1 + 7.2);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(more_dimensions_test, *tolerance(TOL)) {
    Summation sum;
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    shape[2] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 3.14;
    int timestep = 0;
    BOOST_TEST(sum(&tsr, timestep) == 3.14);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(more_dimensions_big_test, *tolerance(TOL)) {
    Summation sum;
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 7.00;
    loc[1] = 1;
    tsr[loc] = 7.01;
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 7.10;
    loc[1] = 1;
    tsr[loc] = 7.11;
    int timestep = 1;
    BOOST_TEST(sum(&tsr, timestep) == 7.00 + 7.01 + 7.10 + 7.11);
    delete[] loc;
}

BOOST_AUTO_TEST_SUITE_END()  // summation_suite
BOOST_AUTO_TEST_SUITE_END()  // fwd_alg_suite

}  // namespace fluoroseq
