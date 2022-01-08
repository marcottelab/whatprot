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
#include "tensor-vector-iterator.h"

// Local project headers:
#include "util/kd-range.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(tensor_suite)
BOOST_AUTO_TEST_SUITE(tensor_vector_iterator_suite)

BOOST_AUTO_TEST_CASE(advance_get_and_end_test) {
    unsigned int order = 3;
    int* stride = new int[order];
    stride[0] = 4;
    stride[1] = 2;
    stride[2] = 1;
    unsigned int vector_dimension = 1;
    KDRange range;
    range.min.resize(order);
    range.max.resize(order);
    range.min[0] = 0;
    range.min[1] = 0;
    range.min[2] = 0;
    range.max[0] = 2;
    range.max[1] = 2;
    range.max[2] = 2;
    unsigned int size = 8;
    double* values = new double[size];
    values[0] = 8000;
    values[1] = 8001;
    values[2] = 8010;
    values[3] = 8011;
    values[4] = 8100;
    values[5] = 8101;
    values[6] = 8110;
    values[7] = 8111;
    TensorVectorIterator itr(
            order, range, range, stride, size, values, vector_dimension);
    Vector* v;
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8000);
    BOOST_TEST((*v)[1] == 8010);
    (*v)[0] = 9000;
    (*v)[1] = 9010;
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8001);
    BOOST_TEST((*v)[1] == 8011);
    (*v)[0] = 9001;
    (*v)[1] = 9011;
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8100);
    BOOST_TEST((*v)[1] == 8110);
    (*v)[0] = 9100;
    (*v)[1] = 9110;
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8101);
    BOOST_TEST((*v)[1] == 8111);
    (*v)[0] = 9101;
    (*v)[1] = 9111;
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == true);
    BOOST_TEST(itr.iter.values[0] == 9000);
    BOOST_TEST(itr.iter.values[1] == 9001);
    BOOST_TEST(itr.iter.values[2] == 9010);
    BOOST_TEST(itr.iter.values[3] == 9011);
    BOOST_TEST(itr.iter.values[4] == 9100);
    BOOST_TEST(itr.iter.values[5] == 9101);
    BOOST_TEST(itr.iter.values[6] == 9110);
    BOOST_TEST(itr.iter.values[7] == 9111);
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_get_and_end_ignore_selected_higher_min_test) {
    unsigned int order = 3;
    int* stride = new int[order];
    stride[0] = 4;
    stride[1] = 2;
    stride[2] = 1;
    unsigned int vector_dimension = 1;
    KDRange itr_range;
    itr_range.min = {1, 1, 1};
    itr_range.max = {2, 2, 2};
    KDRange tsr_range;
    tsr_range.min = {0, 0, 0};
    tsr_range.max = {2, 2, 2};
    unsigned int size = 8;
    double* values = new double[size];
    values[0] = 8000;
    values[1] = 8001;
    values[2] = 8010;
    values[3] = 8011;
    values[4] = 8100;
    values[5] = 8101;
    values[6] = 8110;
    values[7] = 8111;
    TensorVectorIterator itr(order,
                             itr_range,
                             tsr_range,
                             stride,
                             size,
                             values,
                             vector_dimension);
    Vector* v;
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8101);
    BOOST_TEST((*v)[1] == 8111);
    (*v)[0] = 9101;
    (*v)[1] = 9111;
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == true);
    BOOST_TEST(itr.iter.values[0] == 8000);
    BOOST_TEST(itr.iter.values[1] == 8001);
    BOOST_TEST(itr.iter.values[2] == 8010);
    BOOST_TEST(itr.iter.values[3] == 8011);
    BOOST_TEST(itr.iter.values[4] == 8100);
    BOOST_TEST(itr.iter.values[5] == 9101);
    BOOST_TEST(itr.iter.values[6] == 8110);
    BOOST_TEST(itr.iter.values[7] == 9111);
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_SUITE_END()  // tensor_vector_iterator_suite
BOOST_AUTO_TEST_SUITE_END()  // tensor_suite

}  // namespace whatprot
