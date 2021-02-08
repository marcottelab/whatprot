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
#include "vector.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(tensor_suite)
BOOST_AUTO_TEST_SUITE(vector_suite)

BOOST_AUTO_TEST_CASE(constructor_test) {
    int length = 1;
    int stride = 1;
    double* values = new double[length * stride];
    Vector v(length, stride, values);
    BOOST_TEST(v.length == 1);
    BOOST_TEST(v.stride == 1);
    BOOST_TEST(v.values != (void*)NULL);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(bracket_op_test) {
    int length = 1;
    int stride = 1;
    double* values = new double[length];
    values[0] = 42;
    Vector v(length, stride, values);
    BOOST_TEST(v[0] == 42);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(stride_one_test) {
    int length = 5;
    int stride = 1;
    double* values = new double[length * stride];
    values[0] = 420;
    values[1] = 421;
    values[2] = 422;
    values[3] = 423;
    values[4] = 424;
    Vector v(length, stride, values);
    BOOST_TEST(v[0] == 420);
    BOOST_TEST(v[1] == 421);
    BOOST_TEST(v[2] == 422);
    BOOST_TEST(v[3] == 423);
    BOOST_TEST(v[4] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(stride_two_test) {
    int length = 3;
    int stride = 2;
    double* values = new double[length * stride];
    values[0] = 420;
    values[1] = -1;
    values[2] = 422;
    values[3] = -1;
    values[4] = 424;
    values[5] = -1;
    Vector v(length, stride, values);
    BOOST_TEST(v[0] == 420);
    BOOST_TEST(v[1] == 422);
    BOOST_TEST(v[2] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(const_stride_one_test) {
    int length = 5;
    int stride = 1;
    double* values = new double[length * stride];
    values[0] = 420;
    values[1] = 421;
    values[2] = 422;
    values[3] = 423;
    values[4] = 424;
    Vector v(length, stride, values);
    const Vector& cv = v;
    BOOST_TEST(cv[0] == 420);
    BOOST_TEST(cv[1] == 421);
    BOOST_TEST(cv[2] == 422);
    BOOST_TEST(cv[3] == 423);
    BOOST_TEST(cv[4] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(const_stride_two_test) {
    int length = 3;
    int stride = 2;
    double* values = new double[length * stride];
    values[0] = 420;
    values[1] = -1;
    values[2] = 422;
    values[3] = -1;
    values[4] = 424;
    values[5] = -1;
    Vector v(length, stride, values);
    const Vector& cv = v;
    BOOST_TEST(cv[0] == 420);
    BOOST_TEST(cv[1] == 422);
    BOOST_TEST(cv[2] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_SUITE_END()  // vector_suite
BOOST_AUTO_TEST_SUITE_END()  // tensor_suite

}  // namespace whatprot
