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
#include "tensor/const-tensor-iterator.h"
#include "tensor/tensor-iterator.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using std::move;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(tensor_suite)
BOOST_AUTO_TEST_SUITE(tensor_suite)

BOOST_AUTO_TEST_CASE(constructor_order_one_test, *tolerance(TOL)) {
    int order = 1;
    int* shape = new int[order];
    shape[0] = 1;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 1);
    BOOST_TEST(t.strides[0] == 1);
    BOOST_ASSERT(t.size == 1);
    BOOST_TEST(t.values[0] == 0.0);
}

BOOST_AUTO_TEST_CASE(constructor_order_one_bigger_test, *tolerance(TOL)) {
    int order = 1;
    int* shape = new int[order];
    shape[0] = 10;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 10);
    BOOST_TEST(t.strides[0] == 1);
    BOOST_TEST(t.size == 10);
    BOOST_TEST(t.values[0] == 0.0);
    BOOST_TEST(t.values[1] == 0.0);
    BOOST_TEST(t.values[2] == 0.0);
    BOOST_TEST(t.values[3] == 0.0);
    BOOST_TEST(t.values[4] == 0.0);
    BOOST_TEST(t.values[5] == 0.0);
    BOOST_TEST(t.values[6] == 0.0);
    BOOST_TEST(t.values[7] == 0.0);
    BOOST_TEST(t.values[8] == 0.0);
    BOOST_TEST(t.values[9] == 0.0);
}

BOOST_AUTO_TEST_CASE(constructor_order_two_test, *tolerance(TOL)) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 5;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 3);
    BOOST_TEST(t.shape[1] == 5);
    BOOST_TEST(t.strides[0] == 5);
    BOOST_TEST(t.strides[1] == 1);
    BOOST_ASSERT(t.size == 3 * 5);
    BOOST_TEST(t.values[0] == 0.0);
    BOOST_TEST(t.values[1] == 0.0);
    BOOST_TEST(t.values[2] == 0.0);
    BOOST_TEST(t.values[3] == 0.0);
    BOOST_TEST(t.values[4] == 0.0);
    BOOST_TEST(t.values[5] == 0.0);
    BOOST_TEST(t.values[6] == 0.0);
    BOOST_TEST(t.values[7] == 0.0);
    BOOST_TEST(t.values[8] == 0.0);
    BOOST_TEST(t.values[9] == 0.0);
    BOOST_TEST(t.values[10] == 0.0);
    BOOST_TEST(t.values[11] == 0.0);
    BOOST_TEST(t.values[12] == 0.0);
    BOOST_TEST(t.values[13] == 0.0);
    BOOST_TEST(t.values[14] == 0.0);
}

BOOST_AUTO_TEST_CASE(constructor_order_three_test, *tolerance(TOL)) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    shape[2] = 4;
    Tensor t(order, shape);
    delete[] shape;
    BOOST_TEST(t.order == order);
    BOOST_TEST(t.shape[0] == 2);
    BOOST_TEST(t.shape[1] == 3);
    BOOST_TEST(t.shape[2] == 4);
    BOOST_TEST(t.strides[0] == 3 * 4);
    BOOST_TEST(t.strides[1] == 4);
    BOOST_TEST(t.strides[2] == 1);
    BOOST_ASSERT(t.size == 2 * 3 * 4);  // size is 24
    BOOST_TEST(t.values[0] == 0.0);
    BOOST_TEST(t.values[1] == 0.0);
    BOOST_TEST(t.values[2] == 0.0);
    BOOST_TEST(t.values[3] == 0.0);
    BOOST_TEST(t.values[4] == 0.0);
    BOOST_TEST(t.values[5] == 0.0);
    BOOST_TEST(t.values[6] == 0.0);
    BOOST_TEST(t.values[7] == 0.0);
    BOOST_TEST(t.values[8] == 0.0);
    BOOST_TEST(t.values[9] == 0.0);
    BOOST_TEST(t.values[10] == 0.0);
    BOOST_TEST(t.values[11] == 0.0);
    BOOST_TEST(t.values[12] == 0.0);
    BOOST_TEST(t.values[13] == 0.0);
    BOOST_TEST(t.values[14] == 0.0);
    BOOST_TEST(t.values[15] == 0.0);
    BOOST_TEST(t.values[16] == 0.0);
    BOOST_TEST(t.values[17] == 0.0);
    BOOST_TEST(t.values[18] == 0.0);
    BOOST_TEST(t.values[19] == 0.0);
    BOOST_TEST(t.values[20] == 0.0);
    BOOST_TEST(t.values[21] == 0.0);
    BOOST_TEST(t.values[22] == 0.0);
    BOOST_TEST(t.values[23] == 0.0);
}

BOOST_AUTO_TEST_CASE(move_constructor_test, *tolerance(TOL)) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    shape[2] = 4;
    Tensor t1(order, shape);
    delete[] shape;
    Tensor t2(move(t1));
    BOOST_TEST(t1.values == (void*)NULL);
    BOOST_TEST(t1.shape == (void*)NULL);
    BOOST_TEST(t1.strides == (void*)NULL);
    BOOST_TEST(t2.order == order);
    BOOST_TEST(t2.shape[0] == 2);
    BOOST_TEST(t2.shape[1] == 3);
    BOOST_TEST(t2.shape[2] == 4);
    BOOST_TEST(t2.strides[0] == 3 * 4);
    BOOST_TEST(t2.strides[1] == 4);
    BOOST_TEST(t2.strides[2] == 1);
    BOOST_ASSERT(t2.size == 2 * 3 * 4);  // size is 24
    BOOST_TEST(t2.values[0] == 0.0);
    BOOST_TEST(t2.values[1] == 0.0);
    BOOST_TEST(t2.values[2] == 0.0);
    BOOST_TEST(t2.values[3] == 0.0);
    BOOST_TEST(t2.values[4] == 0.0);
    BOOST_TEST(t2.values[5] == 0.0);
    BOOST_TEST(t2.values[6] == 0.0);
    BOOST_TEST(t2.values[7] == 0.0);
    BOOST_TEST(t2.values[8] == 0.0);
    BOOST_TEST(t2.values[9] == 0.0);
    BOOST_TEST(t2.values[10] == 0.0);
    BOOST_TEST(t2.values[11] == 0.0);
    BOOST_TEST(t2.values[12] == 0.0);
    BOOST_TEST(t2.values[13] == 0.0);
    BOOST_TEST(t2.values[14] == 0.0);
    BOOST_TEST(t2.values[15] == 0.0);
    BOOST_TEST(t2.values[16] == 0.0);
    BOOST_TEST(t2.values[17] == 0.0);
    BOOST_TEST(t2.values[18] == 0.0);
    BOOST_TEST(t2.values[19] == 0.0);
    BOOST_TEST(t2.values[20] == 0.0);
    BOOST_TEST(t2.values[21] == 0.0);
    BOOST_TEST(t2.values[22] == 0.0);
    BOOST_TEST(t2.values[23] == 0.0);
}

BOOST_AUTO_TEST_CASE(bracket_op_test, *tolerance(TOL)) {
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

BOOST_AUTO_TEST_CASE(iterator_test, *tolerance(TOL)) {
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
    delete itr;
}

BOOST_AUTO_TEST_CASE(const_iterator_test, *tolerance(TOL)) {
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
    ConstTensorIterator* itr = t.const_iterator();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(itr->get() == 500);
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(itr->get() == 501);
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(itr->get() == 502);
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(itr->get() == 510);
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(itr->get() == 511);
    itr->advance();
    BOOST_TEST(itr->done() == false);
    BOOST_TEST(itr->get() == 512);
    itr->advance();
    BOOST_TEST(itr->done() == true);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 500);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 501);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 502);
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(t[loc] == 510);
    loc[1] = 1;
    BOOST_TEST(t[loc] == 511);
    loc[1] = 2;
    BOOST_TEST(t[loc] == 512);
    delete[] loc;
    delete itr;
}

BOOST_AUTO_TEST_CASE(sum_trivial_test, *tolerance(TOL)) {
    int order = 1;
    int* shape = new int[order];
    shape[0] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    tsr[loc] = 3.14;
    BOOST_TEST(tsr.sum() == 3.14);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(sum_bigger_size_test, *tolerance(TOL)) {
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
    BOOST_TEST(tsr.sum() == 7.0 + 7.1 + 7.2);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(sum_more_dimensions_test, *tolerance(TOL)) {
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
    BOOST_TEST(tsr.sum() == 3.14);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(sum_more_dimensions_big_test, *tolerance(TOL)) {
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
    BOOST_TEST(tsr.sum() == 7.00 + 7.01 + 7.10 + 7.11);
    delete[] loc;
}

BOOST_AUTO_TEST_SUITE_END()  // tensor_suite
BOOST_AUTO_TEST_SUITE_END()  // tensor_suite

}  // namespace whatprot
