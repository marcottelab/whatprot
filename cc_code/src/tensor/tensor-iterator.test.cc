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
#include "tensor-iterator.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(tensor_suite)
BOOST_AUTO_TEST_SUITE(tensor_iterator_suite)

BOOST_AUTO_TEST_CASE(constructor_order_one_test) {
    unsigned int order = 1;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 1;
    min[0] = 0;
    max[0] = 1;
    unsigned int size = 1;
    double* values = new double[size];
    values[0] = 13;
    TensorIterator itr(order, min, max, shape, size, values);
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 1u);
    BOOST_TEST(itr.size == 1u);
    BOOST_TEST(itr.values[0] == 13);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_order_two_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 3;
    min[0] = 0;
    min[1] = 0;
    max[0] = 2;
    max[1] = 3;
    unsigned int size = 6;
    double* values = new double[size];
    values[0] = 600;
    values[1] = 601;
    values[2] = 602;
    values[3] = 610;
    values[4] = 611;
    values[5] = 612;
    TensorIterator itr(order, min, max, shape, size, values);
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 2u);
    BOOST_TEST(itr.shape[1] == 3u);
    BOOST_TEST(itr.size == 6u);
    BOOST_TEST(itr.values[0] == 600);
    BOOST_TEST(itr.values[1] == 601);
    BOOST_TEST(itr.values[2] == 602);
    BOOST_TEST(itr.values[3] == 610);
    BOOST_TEST(itr.values[4] == 611);
    BOOST_TEST(itr.values[5] == 612);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(constructor_order_three_test) {
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    min[0] = 0;
    min[1] = 0;
    min[2] = 0;
    max[0] = 2;
    max[1] = 2;
    max[2] = 2;
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
    TensorIterator itr(order, min, max, shape, size, values);
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 2u);
    BOOST_TEST(itr.shape[1] == 2u);
    BOOST_TEST(itr.shape[2] == 2u);
    BOOST_TEST(itr.size == 8u);
    BOOST_TEST(itr.values[0] == 8000);
    BOOST_TEST(itr.values[1] == 8001);
    BOOST_TEST(itr.values[2] == 8010);
    BOOST_TEST(itr.values[3] == 8011);
    BOOST_TEST(itr.values[4] == 8100);
    BOOST_TEST(itr.values[5] == 8101);
    BOOST_TEST(itr.values[6] == 8110);
    BOOST_TEST(itr.values[7] == 8111);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(reset_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 3;
    min[0] = 0;
    min[1] = 0;
    max[0] = 2;
    max[1] = 3;
    unsigned int size = 6;
    double* values = new double[size];
    TensorIterator itr(order, min, max, shape, size, values);
    itr.loc[0] = 17;
    itr.loc[1] = 19;
    itr.reset();
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 0u);
    BOOST_TEST(itr.index == 0u);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(get_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 3;
    min[0] = 0;
    min[1] = 0;
    max[0] = 2;
    max[1] = 3;
    unsigned int size = 6;
    double* values = new double[size];
    values[0] = 600;
    values[1] = 601;
    values[2] = 602;
    values[3] = 610;
    values[4] = 611;
    values[5] = 612;
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();
    *itr.get() = 42;
    BOOST_TEST(itr.order == order);
    BOOST_TEST(itr.shape[0] == 2u);
    BOOST_TEST(itr.shape[1] == 3u);
    BOOST_TEST(itr.size == 6u);
    BOOST_TEST(itr.values[0] == 42);
    BOOST_TEST(itr.values[1] == 601);
    BOOST_TEST(itr.values[2] == 602);
    BOOST_TEST(itr.values[3] == 610);
    BOOST_TEST(itr.values[4] == 611);
    BOOST_TEST(itr.values[5] == 612);
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 0u);
    BOOST_TEST(itr.index == 0u);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 3;
    min[0] = 0;
    min[1] = 0;
    max[0] = 2;
    max[1] = 3;
    unsigned int size = 6;
    double* values = new double[size];
    values[0] = 500;
    values[1] = 501;
    values[2] = 502;
    values[3] = 510;
    values[4] = 511;
    values[5] = 512;
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();
    BOOST_TEST(*itr.get() == 500);
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 0u);
    BOOST_TEST(itr.index == 0u);
    *itr.get() = 600;
    BOOST_TEST(*itr.get() == 600);
    itr.advance();
    BOOST_TEST(*itr.get() == 501);
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 1u);
    *itr.get() = 601;
    BOOST_TEST(*itr.get() == 601);
    itr.advance();
    BOOST_TEST(*itr.get() == 502);
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 2u);
    *itr.get() = 602;
    BOOST_TEST(*itr.get() == 602);
    itr.advance();
    BOOST_TEST(*itr.get() == 510);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 0u);
    BOOST_TEST(itr.index == 3u);
    *itr.get() = 610;
    BOOST_TEST(*itr.get() == 610);
    itr.advance();
    BOOST_TEST(*itr.get() == 511);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 4u);
    *itr.get() = 611;
    BOOST_TEST(*itr.get() == 611);
    itr.advance();
    BOOST_TEST(*itr.get() == 512);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 5u);
    *itr.get() = 612;
    BOOST_TEST(*itr.get() == 612);
    BOOST_TEST(itr.values[0] == 600);
    BOOST_TEST(itr.values[1] == 601);
    BOOST_TEST(itr.values[2] == 602);
    BOOST_TEST(itr.values[3] == 610);
    BOOST_TEST(itr.values[4] == 611);
    BOOST_TEST(itr.values[5] == 612);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_higher_min_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 4;
    min[0] = 1;
    min[1] = 1;
    max[0] = 3;
    max[1] = 4;
    unsigned int size = 12;
    double* values = new double[size];
    values[0] = 500;
    values[1] = 501;
    values[2] = 502;
    values[3] = 503;
    values[4] = 510;
    values[5] = 511;
    values[6] = 512;
    values[7] = 513;
    values[8] = 520;
    values[9] = 521;
    values[10] = 522;
    values[11] = 523;
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();
    BOOST_TEST(*itr.get() == 511);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 5u);
    *itr.get() = 611;
    BOOST_TEST(*itr.get() == 611);
    itr.advance();
    BOOST_TEST(*itr.get() == 512);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 6u);
    *itr.get() = 612;
    BOOST_TEST(*itr.get() == 612);
    itr.advance();
    BOOST_TEST(*itr.get() == 513);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 3u);
    BOOST_TEST(itr.index == 7u);
    *itr.get() = 613;
    BOOST_TEST(*itr.get() == 613);
    itr.advance();
    BOOST_TEST(*itr.get() == 521);
    BOOST_TEST(itr.loc[0] == 2u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 9u);
    *itr.get() = 621;
    BOOST_TEST(*itr.get() == 621);
    itr.advance();
    BOOST_TEST(*itr.get() == 522);
    BOOST_TEST(itr.loc[0] == 2u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 10u);
    *itr.get() = 622;
    BOOST_TEST(*itr.get() == 622);
    itr.advance();
    BOOST_TEST(*itr.get() == 523);
    BOOST_TEST(itr.loc[0] == 2u);
    BOOST_TEST(itr.loc[1] == 3u);
    BOOST_TEST(itr.index == 11u);
    *itr.get() = 623;
    BOOST_TEST(*itr.get() == 623);
    BOOST_TEST(itr.values[0] == 500);
    BOOST_TEST(itr.values[1] == 501);
    BOOST_TEST(itr.values[2] == 502);
    BOOST_TEST(itr.values[3] == 503);
    BOOST_TEST(itr.values[4] == 510);
    BOOST_TEST(itr.values[5] == 611);
    BOOST_TEST(itr.values[6] == 612);
    BOOST_TEST(itr.values[7] == 613);
    BOOST_TEST(itr.values[8] == 520);
    BOOST_TEST(itr.values[9] == 621);
    BOOST_TEST(itr.values[10] == 622);
    BOOST_TEST(itr.values[11] == 623);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_lower_max_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 4;
    min[0] = 0;
    min[1] = 0;
    max[0] = 2;
    max[1] = 3;
    unsigned int size = 12;
    double* values = new double[size];
    values[0] = 500;
    values[1] = 501;
    values[2] = 502;
    values[3] = 503;
    values[4] = 510;
    values[5] = 511;
    values[6] = 512;
    values[7] = 513;
    values[8] = 520;
    values[9] = 521;
    values[10] = 522;
    values[11] = 523;
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();
    BOOST_TEST(*itr.get() == 500);
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 0u);
    BOOST_TEST(itr.index == 0u);
    *itr.get() = 600;
    BOOST_TEST(*itr.get() == 600);
    itr.advance();
    BOOST_TEST(*itr.get() == 501);
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 1u);
    *itr.get() = 601;
    BOOST_TEST(*itr.get() == 601);
    itr.advance();
    BOOST_TEST(*itr.get() == 502);
    BOOST_TEST(itr.loc[0] == 0u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 2u);
    *itr.get() = 602;
    BOOST_TEST(*itr.get() == 602);
    itr.advance();
    BOOST_TEST(*itr.get() == 510);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 0u);
    BOOST_TEST(itr.index == 4u);
    *itr.get() = 610;
    BOOST_TEST(*itr.get() == 610);
    itr.advance();
    BOOST_TEST(*itr.get() == 511);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 5u);
    *itr.get() = 611;
    BOOST_TEST(*itr.get() == 611);
    itr.advance();
    BOOST_TEST(*itr.get() == 512);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 6u);
    *itr.get() = 612;
    BOOST_TEST(*itr.get() == 612);
    BOOST_TEST(itr.values[0] == 600);
    BOOST_TEST(itr.values[1] == 601);
    BOOST_TEST(itr.values[2] == 602);
    BOOST_TEST(itr.values[3] == 503);
    BOOST_TEST(itr.values[4] == 610);
    BOOST_TEST(itr.values[5] == 611);
    BOOST_TEST(itr.values[6] == 612);
    BOOST_TEST(itr.values[7] == 513);
    BOOST_TEST(itr.values[8] == 520);
    BOOST_TEST(itr.values[9] == 521);
    BOOST_TEST(itr.values[10] == 522);
    BOOST_TEST(itr.values[11] == 523);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(advance_higher_min_lower_max_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 4;
    shape[1] = 5;
    min[0] = 1;
    min[1] = 1;
    max[0] = 3;
    max[1] = 4;
    unsigned int size = 20;
    double* values = new double[size];
    values[0] = 500;
    values[1] = 501;
    values[2] = 502;
    values[3] = 503;
    values[4] = 504;
    values[5] = 510;
    values[6] = 511;
    values[7] = 512;
    values[8] = 513;
    values[9] = 514;
    values[10] = 520;
    values[11] = 521;
    values[12] = 522;
    values[13] = 523;
    values[14] = 524;
    values[15] = 530;
    values[16] = 531;
    values[17] = 532;
    values[18] = 533;
    values[19] = 534;
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();
    BOOST_TEST(*itr.get() == 511);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 6u);
    *itr.get() = 611;
    BOOST_TEST(*itr.get() == 611);
    itr.advance();
    BOOST_TEST(*itr.get() == 512);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 7u);
    *itr.get() = 612;
    BOOST_TEST(*itr.get() == 612);
    itr.advance();
    BOOST_TEST(*itr.get() == 513);
    BOOST_TEST(itr.loc[0] == 1u);
    BOOST_TEST(itr.loc[1] == 3u);
    BOOST_TEST(itr.index == 8u);
    *itr.get() = 613;
    BOOST_TEST(*itr.get() == 613);
    itr.advance();
    BOOST_TEST(*itr.get() == 521);
    BOOST_TEST(itr.loc[0] == 2u);
    BOOST_TEST(itr.loc[1] == 1u);
    BOOST_TEST(itr.index == 11u);
    *itr.get() = 621;
    BOOST_TEST(*itr.get() == 621);
    itr.advance();
    BOOST_TEST(*itr.get() == 522);
    BOOST_TEST(itr.loc[0] == 2u);
    BOOST_TEST(itr.loc[1] == 2u);
    BOOST_TEST(itr.index == 12u);
    *itr.get() = 622;
    BOOST_TEST(*itr.get() == 622);
    itr.advance();
    BOOST_TEST(*itr.get() == 523);
    BOOST_TEST(itr.loc[0] == 2u);
    BOOST_TEST(itr.loc[1] == 3u);
    BOOST_TEST(itr.index == 13u);
    *itr.get() = 623;
    BOOST_TEST(*itr.get() == 623);
    BOOST_TEST(itr.values[0] == 500);
    BOOST_TEST(itr.values[1] == 501);
    BOOST_TEST(itr.values[2] == 502);
    BOOST_TEST(itr.values[3] == 503);
    BOOST_TEST(itr.values[4] == 504);
    BOOST_TEST(itr.values[5] == 510);
    BOOST_TEST(itr.values[6] == 611);
    BOOST_TEST(itr.values[7] == 612);
    BOOST_TEST(itr.values[8] == 613);
    BOOST_TEST(itr.values[9] == 514);
    BOOST_TEST(itr.values[10] == 520);
    BOOST_TEST(itr.values[11] == 621);
    BOOST_TEST(itr.values[12] == 622);
    BOOST_TEST(itr.values[13] == 623);
    BOOST_TEST(itr.values[14] == 524);
    BOOST_TEST(itr.values[15] == 530);
    BOOST_TEST(itr.values[16] == 531);
    BOOST_TEST(itr.values[17] == 532);
    BOOST_TEST(itr.values[18] == 533);
    BOOST_TEST(itr.values[19] == 534);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(done_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 3;
    min[0] = 0;
    min[1] = 0;
    max[0] = 2;
    max[1] = 3;
    unsigned int size = 6;
    double* values = new double[size];
    TensorIterator itr(order, min, max, shape, size, values);
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
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(done_higher_min_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 4;
    min[0] = 1;
    min[1] = 1;
    max[0] = 3;
    max[1] = 4;
    unsigned int size = 12;
    double* values = new double[size];
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();  // at (1, 1)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (1, 2)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (1, 3)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (2, 1)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (2, 2)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (2, 3)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // should be done
    BOOST_TEST(itr.done() == true);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(done_lower_max_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 4;
    min[0] = 0;
    min[1] = 0;
    max[0] = 2;
    max[1] = 3;
    unsigned int size = 12;
    double* values = new double[size];
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();  // at (0, 0)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (0, 1)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (0, 2)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (1, 0)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (1, 1)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (1, 2)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // should be done
    BOOST_TEST(itr.done() == true);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(done_higher_min_lower_max_test) {
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 4;
    shape[1] = 5;
    min[0] = 1;
    min[1] = 1;
    max[0] = 3;
    max[1] = 4;
    unsigned int size = 20;
    double* values = new double[size];
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();  // at (1, 1)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (1, 2)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (1, 3)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (2, 1)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (2, 2)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // at (2, 3)
    BOOST_TEST(itr.done() == false);
    itr.advance();  // should be done
    BOOST_TEST(itr.done() == true);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_CASE(order_three_size_one_test) {
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    unsigned int* min = new unsigned int[order];
    unsigned int* max = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    shape[2] = 1;
    min[0] = 0;
    min[1] = 0;
    min[2] = 0;
    max[0] = 1;
    max[1] = 1;
    max[2] = 1;
    unsigned int size = 1;
    double* values = new double[size];
    values[0] = 13;
    TensorIterator itr(order, min, max, shape, size, values);
    itr.reset();
    BOOST_TEST(itr.done() == false);
    BOOST_TEST(*itr.get() == 13);
    itr.advance();
    BOOST_TEST(itr.done() == true);
    delete[] shape;
    delete[] min;
    delete[] max;
    delete[] values;
}

BOOST_AUTO_TEST_SUITE_END()  // tensor_iterator_suite
BOOST_AUTO_TEST_SUITE_END()  // tensor_suite

}  // namespace whatprot
