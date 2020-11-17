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
#include "tensor.h"

// Standard C++ library headers:
#include <utility>

// Local project headers:
#include "tensor/tensor_iterator.h"

namespace fluoroseq {

namespace {
using std::move;
}  // namespace

BOOST_AUTO_TEST_SUITE(tensor_suite);
BOOST_AUTO_TEST_SUITE(tensor_suite);

BOOST_AUTO_TEST_CASE(constructor_order_one_test) {
    int order = 1;
    int* shape = new int[order];
    shape[0] = 1;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.tensor_iterator != (void*)NULL);
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 1);
    BOOST_TEST(t.strides[0] == 1);
    BOOST_TEST(t.size == 1);
}

BOOST_AUTO_TEST_CASE(constructor_order_one_bigger_test) {
    int order = 1;
    int* shape = new int[order];
    shape[0] = 10;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.tensor_iterator != (void*)NULL);
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 10);
    BOOST_TEST(t.strides[0] == 1);
    BOOST_TEST(t.size == 10);
}

BOOST_AUTO_TEST_CASE(constructor_order_two_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 5;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.tensor_iterator != (void*)NULL);
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 3);
    BOOST_TEST(t.shape[1] == 5);
    BOOST_TEST(t.strides[0] == 5);
    BOOST_TEST(t.strides[1] == 1);
    BOOST_TEST(t.size == 3 * 5);
}

BOOST_AUTO_TEST_CASE(constructor_order_three_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 5;
    shape[2] = 7;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.tensor_iterator != (void*)NULL);
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 3);
    BOOST_TEST(t.shape[1] == 5);
    BOOST_TEST(t.shape[2] == 7);
    BOOST_TEST(t.strides[0] == 5 * 7);
    BOOST_TEST(t.strides[1] == 7);
    BOOST_TEST(t.strides[2] == 1);
    BOOST_TEST(t.size == 3 * 5 * 7);
}

BOOST_AUTO_TEST_CASE(move_constructor_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 5;
    shape[2] = 7;
    Tensor t1(order, shape);
    delete[] shape;
    Tensor t2(move(t1));
    BOOST_TEST(t1.tensor_iterator == (void*)NULL);
    BOOST_TEST(t1.values == (void*)NULL);
    BOOST_TEST(t1.shape == (void*)NULL);
    BOOST_TEST(t1.strides == (void*)NULL);
    BOOST_TEST(t2.order == order);
    BOOST_TEST(t2.shape[0] == 3);
    BOOST_TEST(t2.shape[1] == 5);
    BOOST_TEST(t2.shape[2] == 7);
    BOOST_TEST(t2.strides[0] == 5 * 7);
    BOOST_TEST(t2.strides[1] == 7);
    BOOST_TEST(t2.strides[2] == 1);
    BOOST_TEST(t2.size == 3 * 5 * 7);
}

BOOST_AUTO_TEST_CASE(bracket_op_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    Tensor t(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    t[loc] = 600;
    loc[1] = 1;
    t[loc] = 601;
    loc[1] = 2;
    t[loc] = 602;
    loc[0] = 1;
    loc[1] = 0;
    t[loc] = 610;
    loc[1] = 1;
    t[loc] = 611;
    loc[1] = 2;
    t[loc] = 612;
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 600);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 601);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 602);
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 610);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 611);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 612);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(iterator_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    Tensor t(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    t[loc] = 500;
    loc[1] = 1;
    t[loc] = 501;
    loc[1] = 2;
    t[loc] = 502;
    loc[0] = 1;
    loc[1] = 0;
    t[loc] = 510;
    loc[1] = 1;
    t[loc] = 511;
    loc[1] = 2;
    t[loc] = 512;
    TensorIterator* itr = t.iterator();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 500);
    *itr->get() = 600;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 501);
    *itr->get() = 601;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 502);
    *itr->get() = 602;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 510);
    *itr->get() = 610;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 511);
    *itr->get() = 611;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 512);
    *itr->get() = 612;
    itr->advance();
    BOOST_TEST(itr->done() == true);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 600);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 601);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 602);
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 610);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 611);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 612);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(iterator_reuse_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    Tensor t(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    t[loc] = 500;
    loc[1] = 1;
    t[loc] = 501;
    loc[1] = 2;
    t[loc] = 502;
    loc[0] = 1;
    loc[1] = 0;
    t[loc] = 510;
    loc[1] = 1;
    t[loc] = 511;
    loc[1] = 2;
    t[loc] = 512;
    TensorIterator* itr;
    itr = t.iterator();
    itr->advance();
    itr->advance();
    itr->advance();
    itr->advance();
    itr = t.iterator();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 500);
    *itr->get() = 600;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 501);
    *itr->get() = 601;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 502);
    *itr->get() = 602;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 510);
    *itr->get() = 610;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 511);
    *itr->get() = 611;
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(*itr->get() == 512);
    *itr->get() = 612;
    itr->advance();
    BOOST_TEST(itr->done() == true);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 600);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 601);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 602);
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 610);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 611);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 612);
    delete[] loc;
}

BOOST_AUTO_TEST_SUITE_END();  // tensor_suite
BOOST_AUTO_TEST_SUITE_END();  // tensor_suite

}  // namespace fluoroseq
