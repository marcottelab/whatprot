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
#include "const-tensor-vector-iterator.h"

// Local project headers:
#include "util/kd-range.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(tensor_suite)
BOOST_AUTO_TEST_SUITE(const_tensor_vector_iterator_suite)

BOOST_AUTO_TEST_CASE(constructor_order_one_test) {
    unsigned int order = 1;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    unsigned int* stride = new unsigned int[order];
    stride[0] = 1;
    unsigned int vector_dimension = 0;
    KDRange range;
    range.min.resize(order);
    range.max.resize(order);
    range.min[0] = 0;
    range.max[0] = 1;
    unsigned int size = 1;
    double* values = new double[size];
    values[0] = 13;
    ConstTensorVectorIterator itr(order, range, shape, stride, size, values, vector_dimension);
    BOOST_TEST(itr.iter.order == order);
    BOOST_TEST(itr.iter.shape[0] == 1u);
    BOOST_TEST(itr.iter.size == 1u);
    BOOST_TEST(itr.iter.values[0] == 13);
    BOOST_TEST(itr.vector_length == 1u);
    BOOST_TEST(itr.vector_stride == 1u);
    BOOST_TEST(itr.modified_range.min[0] == 0u);
    BOOST_TEST(itr.iter.range.min[0] == 0u);
    BOOST_TEST(itr.iter.loc[0] == 0u);
    BOOST_TEST(itr.modified_range.max[0] == 1u);
    BOOST_TEST(itr.iter.range.max[0] == 1u);
    delete[] shape;
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_modifies_range_test) {
    unsigned int order = 1;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    unsigned int* stride = new unsigned int[order];
    stride[0] = 1;
    unsigned int vector_dimension = 0;
    KDRange range;
    range.min.resize(order);
    range.max.resize(order);
    range.min[0] = 2;
    range.max[0] = 3;
    unsigned int size = 3;
    double* values = new double[size];
    values[0] = 13;
    values[1] = 14;
    values[2] = 15;
    ConstTensorVectorIterator itr(order, range, shape, stride, size, values, vector_dimension);
    BOOST_TEST(itr.iter.order == order);
    BOOST_TEST(itr.iter.shape[0] == 3u);
    BOOST_TEST(itr.iter.size == 3u);
    BOOST_TEST(itr.iter.values[0] == 13);
    BOOST_TEST(itr.iter.values[1] == 14);
    BOOST_TEST(itr.iter.values[2] == 15);
    BOOST_TEST(itr.vector_length == 3u);
    BOOST_TEST(itr.vector_stride == 1u);
    BOOST_TEST(itr.modified_range.min[0] == 0u);
    BOOST_TEST(itr.iter.range.min[0] == 0u);
    BOOST_TEST(itr.iter.loc[0] == 0u);
    BOOST_TEST(itr.modified_range.max[0] == 1u);
    BOOST_TEST(itr.iter.range.max[0] == 1u);
    delete[] shape;
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_order_two_vectorize_dim_0_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 3;
    unsigned int* stride = new unsigned int[order];
    stride[0] = 3;
    stride[1] = 1;
    unsigned int vector_dimension = 0;
    KDRange range;
    range.min.resize(order);
    range.max.resize(order);
    range.min[0] = 1;
    range.min[1] = 0;
    range.max[0] = 2;
    range.max[1] = 3;
    unsigned int size = 6;
    double* values = new double[size];
    values[0] = 600;
    values[1] = 601;
    values[2] = 602;
    values[3] = 610;
    values[4] = 611;
    values[5] = 612;
    ConstTensorVectorIterator itr(order, range, shape, stride, size, values, vector_dimension);
    BOOST_TEST(itr.iter.order == order);
    BOOST_TEST(itr.iter.shape[0] == 2u);
    BOOST_TEST(itr.iter.shape[1] == 3u);
    BOOST_TEST(itr.iter.size == 6u);
    BOOST_TEST(itr.iter.values[0] == 600);
    BOOST_TEST(itr.iter.values[1] == 601);
    BOOST_TEST(itr.iter.values[2] == 602);
    BOOST_TEST(itr.iter.values[3] == 610);
    BOOST_TEST(itr.iter.values[4] == 611);
    BOOST_TEST(itr.iter.values[5] == 612);
    BOOST_TEST(itr.vector_length == 2u);
    BOOST_TEST(itr.vector_stride == 3u);
    BOOST_TEST(itr.modified_range.min[0] == 0u);
    BOOST_TEST(itr.iter.range.min[0] == 0u);
    BOOST_TEST(itr.iter.loc[0] == 0u);
    BOOST_TEST(itr.modified_range.max[0] == 1u);
    BOOST_TEST(itr.iter.range.max[0] == 1u);
    BOOST_TEST(itr.modified_range.min[1] == 0u);
    BOOST_TEST(itr.iter.range.min[1] == 0u);
    BOOST_TEST(itr.iter.loc[1] == 0u);
    BOOST_TEST(itr.modified_range.max[1] == 3u);
    BOOST_TEST(itr.iter.range.max[1] == 3u);
    delete[] shape;
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_order_two_vectorize_dim_1_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 3;
    unsigned int* stride = new unsigned int[order];
    stride[0] = 3;
    stride[1] = 1;
    unsigned int vector_dimension = 1;
    KDRange range;
    range.min.resize(order);
    range.max.resize(order);
    range.min[0] = 0;
    range.min[1] = 2;
    range.max[0] = 2;
    range.max[1] = 3;
    unsigned int size = 6;
    double* values = new double[size];
    values[0] = 600;
    values[1] = 601;
    values[2] = 602;
    values[3] = 610;
    values[4] = 611;
    values[5] = 612;
    ConstTensorVectorIterator itr(order, range, shape, stride, size, values, vector_dimension);
    BOOST_TEST(itr.iter.order == order);
    BOOST_TEST(itr.iter.shape[0] == 2u);
    BOOST_TEST(itr.iter.shape[1] == 3u);
    BOOST_TEST(itr.iter.size == 6u);
    BOOST_TEST(itr.iter.values[0] == 600);
    BOOST_TEST(itr.iter.values[1] == 601);
    BOOST_TEST(itr.iter.values[2] == 602);
    BOOST_TEST(itr.iter.values[3] == 610);
    BOOST_TEST(itr.iter.values[4] == 611);
    BOOST_TEST(itr.iter.values[5] == 612);
    BOOST_TEST(itr.vector_length == 3u);
    BOOST_TEST(itr.vector_stride == 1u);
    BOOST_TEST(itr.modified_range.min[0] == 0u);
    BOOST_TEST(itr.iter.range.min[0] == 0u);
    BOOST_TEST(itr.iter.loc[0] == 0u);
    BOOST_TEST(itr.modified_range.max[0] == 2u);
    BOOST_TEST(itr.iter.range.max[0] == 2u);
    BOOST_TEST(itr.modified_range.min[1] == 0u);
    BOOST_TEST(itr.iter.range.min[1] == 0u);
    BOOST_TEST(itr.iter.loc[1] == 0u);
    BOOST_TEST(itr.modified_range.max[1] == 1u);
    BOOST_TEST(itr.iter.range.max[1] == 1u);
    delete[] shape;
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_order_three_test) {
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    unsigned int* stride = new unsigned int[order];
    stride[0] = 4;
    stride[1] = 2;
    stride[2] = 1;
    unsigned int vector_dimension = 1;
    KDRange range;
    range.min.resize(order);
    range.max.resize(order);
    range.min[0] = 0;
    range.min[1] = 1;
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
    ConstTensorVectorIterator itr(order, range, shape, stride, size, values, vector_dimension);
    BOOST_TEST(itr.iter.order == order);
    BOOST_TEST(itr.iter.shape[0] == 2u);
    BOOST_TEST(itr.iter.shape[1] == 2u);
    BOOST_TEST(itr.iter.shape[2] == 2u);
    BOOST_TEST(itr.iter.size == 8u);
    BOOST_TEST(itr.iter.values[0] == 8000);
    BOOST_TEST(itr.iter.values[1] == 8001);
    BOOST_TEST(itr.iter.values[2] == 8010);
    BOOST_TEST(itr.iter.values[3] == 8011);
    BOOST_TEST(itr.iter.values[4] == 8100);
    BOOST_TEST(itr.iter.values[5] == 8101);
    BOOST_TEST(itr.iter.values[6] == 8110);
    BOOST_TEST(itr.iter.values[7] == 8111);
    BOOST_TEST(itr.vector_length == 2u);
    BOOST_TEST(itr.vector_stride == 2u);
    BOOST_TEST(itr.modified_range.min[0] == 0u);
    BOOST_TEST(itr.iter.range.min[0] == 0u);
    BOOST_TEST(itr.iter.loc[0] == 0u);
    BOOST_TEST(itr.modified_range.max[0] == 2u);
    BOOST_TEST(itr.iter.range.max[0] == 2u);
    BOOST_TEST(itr.modified_range.min[1] == 0u);
    BOOST_TEST(itr.iter.range.min[1] == 0u);
    BOOST_TEST(itr.iter.loc[1] == 0u);
    BOOST_TEST(itr.modified_range.max[1] == 1u);
    BOOST_TEST(itr.iter.range.max[1] == 1u);
    BOOST_TEST(itr.modified_range.min[2] == 0u);
    BOOST_TEST(itr.iter.range.min[2] == 0u);
    BOOST_TEST(itr.iter.loc[2] == 0u);
    BOOST_TEST(itr.modified_range.max[2] == 2u);
    BOOST_TEST(itr.iter.range.max[2] == 2u);
    delete[] shape;
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_get_and_end_test) {
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    unsigned int* stride = new unsigned int[order];
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
    ConstTensorVectorIterator itr(order, range, shape, stride, size, values, vector_dimension);
    const Vector* v;
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8000);
    BOOST_TEST((*v)[1] == 8010);
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8001);
    BOOST_TEST((*v)[1] == 8011);
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8100);
    BOOST_TEST((*v)[1] == 8110);
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8101);
    BOOST_TEST((*v)[1] == 8111);
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == true);
    BOOST_TEST(itr.iter.values[0] == 8000);
    BOOST_TEST(itr.iter.values[1] == 8001);
    BOOST_TEST(itr.iter.values[2] == 8010);
    BOOST_TEST(itr.iter.values[3] == 8011);
    BOOST_TEST(itr.iter.values[4] == 8100);
    BOOST_TEST(itr.iter.values[5] == 8101);
    BOOST_TEST(itr.iter.values[6] == 8110);
    BOOST_TEST(itr.iter.values[7] == 8111);
    delete[] shape;
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_get_and_end_ignore_selected_higher_min_test) {
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    unsigned int* stride = new unsigned int[order];
    stride[0] = 4;
    stride[1] = 2;
    stride[2] = 1;
    unsigned int vector_dimension = 1;
    KDRange range;
    range.min.resize(order);
    range.max.resize(order);
    range.min[0] = 1;
    range.min[1] = 1;
    range.min[2] = 1;
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
    ConstTensorVectorIterator itr(order, range, shape, stride, size, values, vector_dimension);
    const Vector* v;
    BOOST_TEST(itr.done() == false);
    v = itr.get();
    BOOST_TEST(v->length == 2u);
    BOOST_TEST(v->stride == 2);
    BOOST_TEST((*v)[0] == 8101);
    BOOST_TEST((*v)[1] == 8111);
    delete v;
    itr.advance();
    BOOST_TEST(itr.done() == true);
    BOOST_TEST(itr.iter.values[0] == 8000);
    BOOST_TEST(itr.iter.values[1] == 8001);
    BOOST_TEST(itr.iter.values[2] == 8010);
    BOOST_TEST(itr.iter.values[3] == 8011);
    BOOST_TEST(itr.iter.values[4] == 8100);
    BOOST_TEST(itr.iter.values[5] == 8101);
    BOOST_TEST(itr.iter.values[6] == 8110);
    BOOST_TEST(itr.iter.values[7] == 8111);
    delete[] shape;
    delete[] stride;
    delete[] values;
}

BOOST_AUTO_TEST_SUITE_END()  // const_tensor_vector_iterator_suite
BOOST_AUTO_TEST_SUITE_END()  // tensor_suite

}  // namespace whatprot
