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
    unsigned int v_min = 0;
    unsigned int v_max = 1;
    int stride = 1;
    double* values = new double[v_max * stride];
    Vector v(v_min, v_max, stride, values);
    BOOST_TEST(v.min == 0u);
    BOOST_TEST(v.max == 1u);
    BOOST_TEST(v.stride == 1);
    BOOST_TEST(v.values != (void*)NULL);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(bracket_op_test) {
    unsigned int v_min = 0;
    unsigned int v_max = 1;
    int stride = 1;
    double* values = new double[v_max];
    values[0] = 42;
    Vector v(v_min, v_max, stride, values);
    BOOST_TEST(v[0] == 42);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(stride_one_test) {
    unsigned int v_min = 0;
    unsigned int v_max = 5;
    int stride = 1;
    double* values = new double[v_max * stride];
    values[0] = 420;
    values[1] = 421;
    values[2] = 422;
    values[3] = 423;
    values[4] = 424;
    Vector v(v_min, v_max, stride, values);
    BOOST_TEST(v[0] == 420);
    BOOST_TEST(v[1] == 421);
    BOOST_TEST(v[2] == 422);
    BOOST_TEST(v[3] == 423);
    BOOST_TEST(v[4] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(offset_tensor_test) {
    unsigned int v_min = 1;
    unsigned int v_max = 5;
    int stride = 1;
    double* values = new double[(v_max - v_min) * stride];
    values[0] = 421;
    values[1] = 422;
    values[2] = 423;
    values[3] = 424;
    Vector v(v_min, v_max, stride, values);
    BOOST_TEST(v[1] == 421);
    BOOST_TEST(v[2] == 422);
    BOOST_TEST(v[3] == 423);
    BOOST_TEST(v[4] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(stride_two_test) {
    unsigned int v_min = 0;
    unsigned int v_max = 3;
    int stride = 2;
    double* values = new double[v_max * stride];
    values[0] = 420;
    values[1] = -1;
    values[2] = 422;
    values[3] = -1;
    values[4] = 424;
    values[5] = -1;
    Vector v(v_min, v_max, stride, values);
    BOOST_TEST(v[0] == 420);
    BOOST_TEST(v[1] == 422);
    BOOST_TEST(v[2] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(mutable_values_test) {
    unsigned int v_min = 0;
    unsigned int v_max = 5;
    int stride = 1;
    double* values = new double[v_max * stride];
    values[0] = 420;
    values[1] = 421;
    values[2] = 422;
    values[3] = 423;
    values[4] = 424;
    Vector v(v_min, v_max, stride, values);
    v[0] = 720;
    v[1] = 721;
    v[2] = 722;
    v[3] = 723;
    v[4] = 724;
    BOOST_TEST(v[0] == 720);
    BOOST_TEST(v[1] == 721);
    BOOST_TEST(v[2] == 722);
    BOOST_TEST(v[3] == 723);
    BOOST_TEST(v[4] == 724);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(mutable_values_offset_tensor_test) {
    unsigned int v_min = 1;
    unsigned int v_max = 5;
    int stride = 1;
    double* values = new double[(v_max - v_min) * stride];
    values[0] = 421;
    values[1] = 422;
    values[2] = 423;
    values[3] = 424;
    Vector v(v_min, v_max, stride, values);
    v[1] = 721;
    v[2] = 722;
    v[3] = 723;
    v[4] = 724;
    BOOST_TEST(v[1] == 721);
    BOOST_TEST(v[2] == 722);
    BOOST_TEST(v[3] == 723);
    BOOST_TEST(v[4] == 724);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(const_stride_one_test) {
    unsigned int v_min = 0;
    unsigned int v_max = 5;
    int stride = 1;
    double* values = new double[v_max * stride];
    values[0] = 420;
    values[1] = 421;
    values[2] = 422;
    values[3] = 423;
    values[4] = 424;
    Vector v(v_min, v_max, stride, values);
    const Vector& cv = v;
    BOOST_TEST(cv[0] == 420);
    BOOST_TEST(cv[1] == 421);
    BOOST_TEST(cv[2] == 422);
    BOOST_TEST(cv[3] == 423);
    BOOST_TEST(cv[4] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_CASE(const_stride_two_test) {
    unsigned int v_min = 0;
    unsigned int v_max = 3;
    int stride = 2;
    double* values = new double[v_max * stride];
    values[0] = 420;
    values[1] = -1;
    values[2] = 422;
    values[3] = -1;
    values[4] = 424;
    values[5] = -1;
    Vector v(v_min, v_max, stride, values);
    const Vector& cv = v;
    BOOST_TEST(cv[0] == 420);
    BOOST_TEST(cv[1] == 422);
    BOOST_TEST(cv[2] == 424);
    delete[] values;
}

BOOST_AUTO_TEST_SUITE_END()  // vector_suite
BOOST_AUTO_TEST_SUITE_END()  // tensor_suite

}  // namespace whatprot
