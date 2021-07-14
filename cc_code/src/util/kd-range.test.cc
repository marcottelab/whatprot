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
#include "kd-range.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(util_suite)
BOOST_AUTO_TEST_SUITE(kd_box_range_suite)

BOOST_AUTO_TEST_CASE(intersect_test) {
    KDRange x;
    x.min.resize(3);
    x.min[0] = 2;
    x.min[1] = 5;
    x.min[2] = 7;
    x.max.resize(3);
    x.max[0] = 4;
    x.max[1] = 9;
    x.max[2] = 11;
    KDRange y;
    y.min.resize(3);
    y.min[0] = 1;
    y.min[1] = 7;
    y.min[2] = 6;
    y.max.resize(3);
    y.max[0] = 5;
    y.max[1] = 13;
    y.max[2] = 10;
    KDRange z = x.intersect(y);
    BOOST_TEST(z.min[0] == 2u);
    BOOST_TEST(z.min[1] == 7u);
    BOOST_TEST(z.min[2] == 7u);
    BOOST_TEST(z.max[0] == 4u);
    BOOST_TEST(z.max[1] == 9u);
    BOOST_TEST(z.max[2] == 10u);
}

BOOST_AUTO_TEST_CASE(is_empty_false_test) {
    KDRange x;
    x.min.resize(3);
    x.min[0] = 2;
    x.min[1] = 5;
    x.min[2] = 7;
    x.max.resize(3);
    x.max[0] = 4;
    x.max[1] = 8;
    x.max[2] = 11;
    BOOST_TEST(x.is_empty() == false);
}

BOOST_AUTO_TEST_CASE(is_empty_true_test) {
    KDRange x;
    x.min.resize(3);
    x.min[0] = 2;
    x.min[1] = 5;
    x.min[2] = 7;
    x.max.resize(3);
    x.max[0] = 4;
    x.max[1] = 5;
    x.max[2] = 11;
    BOOST_TEST(x.is_empty() == true);
}

BOOST_AUTO_TEST_SUITE_END()  // kd_box_range_suite
BOOST_AUTO_TEST_SUITE_END()  // util_suite

}  // namespace whatprot
