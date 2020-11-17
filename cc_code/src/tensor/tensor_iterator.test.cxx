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
#include "tensor_iterator.h"

namespace fluoroseq {

BOOST_AUTO_TEST_SUITE(tensor_suite);
BOOST_AUTO_TEST_SUITE(tensor_iterator_suite);

BOOST_AUTO_TEST_CASE(constructor_order_one_test) {
    int order = 1;
    int* shape = new int[order];
    shape[0] = 1;
    int size = 1;
    double* values = new double[size];
    values[0] = 13;
    TensorIterator itr(order, shape, size, values);
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 1);
    BOOST_TEST(itr.size == 1);
    BOOST_TEST(itr.values[0] == 13);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_order_two_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    int size = 6;
    double* values = new double[size];
    values[0] = 600;
    values[1] = 601;
    values[2] = 602;
    values[3] = 610;
    values[4] = 611;
    values[5] = 612;
    TensorIterator itr(order, shape, size, values);
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 2);
    BOOST_TEST(itr.shape[1] == 3);
    BOOST_TEST(itr.size == 6);
    BOOST_TEST(itr.values[0] == 600);
    BOOST_TEST(itr.values[1] == 601);
    BOOST_TEST(itr.values[2] == 602);
    BOOST_TEST(itr.values[3] == 610);
    BOOST_TEST(itr.values[4] == 611);
    BOOST_TEST(itr.values[5] == 612);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_order_three_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    int size = 8;
    double* values = new double[size];
    values[0] = 8000;
    values[1] = 8001;
    values[2] = 8010;
    values[3] = 8011;
    values[4] = 8100;
    values[5] = 8101;
    values[6] = 8110;
    values[7] = 8111;
    TensorIterator itr(order, shape, size, values);
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 2);
    BOOST_TEST(itr.shape[1] == 2);
    BOOST_TEST(itr.shape[2] == 2);
    BOOST_TEST(itr.size == 8);
    BOOST_TEST(itr.values[0] == 8000);
    BOOST_TEST(itr.values[1] == 8001);
    BOOST_TEST(itr.values[2] == 8010);
    BOOST_TEST(itr.values[3] == 8011);
    BOOST_TEST(itr.values[4] == 8100);
    BOOST_TEST(itr.values[5] == 8101);
    BOOST_TEST(itr.values[6] == 8110);
    BOOST_TEST(itr.values[7] == 8111);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(reset_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    int size = 6;
    double* values = new double[size];
    TensorIterator itr(order, shape, size, values);
    itr.loc[0] = 17;
    itr.loc[1] = 19;
    itr.reset();
    BOOST_TEST(itr.loc[0] == 0);
    BOOST_TEST(itr.loc[1] == 0);
    BOOST_TEST(itr.index == 0);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(get_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    int size = 6;
    double* values = new double[size];
    values[0] = 600;
    values[1] = 601;
    values[2] = 602;
    values[3] = 610;
    values[4] = 611;
    values[5] = 612;
    TensorIterator itr(order, shape, size, values);
    itr.reset();
    *itr.get() = 42;
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 2);
    BOOST_TEST(itr.shape[1] == 3);
    BOOST_TEST(itr.size == 6);
    BOOST_TEST(itr.values[0] == 42);
    BOOST_TEST(itr.values[1] == 601);
    BOOST_TEST(itr.values[2] == 602);
    BOOST_TEST(itr.values[3] == 610);
    BOOST_TEST(itr.values[4] == 611);
    BOOST_TEST(itr.values[5] == 612);
    BOOST_TEST(itr.loc[0] == 0);
    BOOST_TEST(itr.loc[1] == 0);
    BOOST_TEST(itr.index == 0);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    int size = 6;
    double* values = new double[size];
    values[0] = 500;
    values[1] = 501;
    values[2] = 502;
    values[3] = 510;
    values[4] = 511;
    values[5] = 512;
    TensorIterator itr(order, shape, size, values);
    itr.reset();
    BOOST_TEST(*itr.get() == 500);
    BOOST_TEST(itr.loc[0] == 0);
    BOOST_TEST(itr.loc[1] == 0);
    BOOST_TEST(itr.index == 0);
    *itr.get() = 600;
    BOOST_TEST(*itr.get() == 600);
    itr.advance();
    BOOST_TEST(*itr.get() == 501);
    BOOST_TEST(itr.loc[0] == 0);
    BOOST_TEST(itr.loc[1] == 1);
    BOOST_TEST(itr.index == 1);
    *itr.get() = 601;
    BOOST_TEST(*itr.get() == 601);
    itr.advance();
    BOOST_TEST(*itr.get() == 502);
    BOOST_TEST(itr.loc[0] == 0);
    BOOST_TEST(itr.loc[1] == 2);
    BOOST_TEST(itr.index == 2);
    *itr.get() = 602;
    BOOST_TEST(*itr.get() == 602);
    itr.advance();
    BOOST_TEST(*itr.get() == 510);
    BOOST_TEST(itr.loc[0] == 1);
    BOOST_TEST(itr.loc[1] == 0);
    BOOST_TEST(itr.index == 3);
    *itr.get() = 610;
    BOOST_TEST(*itr.get() == 610);
    itr.advance();
    BOOST_TEST(*itr.get() == 511);
    BOOST_TEST(itr.loc[0] == 1);
    BOOST_TEST(itr.loc[1] == 1);
    BOOST_TEST(itr.index == 4);
    *itr.get() = 611;
    BOOST_TEST(*itr.get() == 611);
    itr.advance();
    BOOST_TEST(*itr.get() == 512);
    BOOST_TEST(itr.loc[0] == 1);
    BOOST_TEST(itr.loc[1] == 2);
    BOOST_TEST(itr.index == 5);
    *itr.get() = 612;
    BOOST_TEST(*itr.get() == 612);
    BOOST_TEST(itr.values[0] == 600);
    BOOST_TEST(itr.values[1] == 601);
    BOOST_TEST(itr.values[2] == 602);
    BOOST_TEST(itr.values[3] == 610);
    BOOST_TEST(itr.values[4] == 611);
    BOOST_TEST(itr.values[5] == 612);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(done_test) {
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 3;
    int size = 6;
    double* values = new double[size];
    TensorIterator itr(order, shape, size, values);
    itr.reset();
    BOOST_TEST(itr.done() == false);
    itr.advance();
    BOOST_TEST(itr.done() == false);
    itr.advance();
    BOOST_TEST(itr.done() == false);
    itr.advance();
    BOOST_TEST(itr.done() == false);
    itr.advance();
    BOOST_TEST(itr.done() == false);
    itr.advance();
    BOOST_TEST(itr.done() == false);
    itr.advance();
    BOOST_TEST(itr.done() == true);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(order_three_size_one_test) {
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    shape[2] = 1;
    int size = 1;
    double* values = new double[size];
    values[0] = 13;
    TensorIterator itr(order, shape, size, values);
    itr.reset();
    BOOST_TEST(itr.done() == false);
    BOOST_TEST(*itr.get() == 13);
    itr.advance();
    BOOST_TEST(itr.done() == true);
    delete[] shape;
    delete[] values;
}

BOOST_AUTO_TEST_SUITE_END();  // tensor_iterator_suite
BOOST_AUTO_TEST_SUITE_END();  // tensor_suite

}  // namespace fluoroseq
