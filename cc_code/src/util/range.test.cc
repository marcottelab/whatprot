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
#include "util/range.h"

namespace fluoroseq {

BOOST_AUTO_TEST_SUITE(util_suite)
BOOST_AUTO_TEST_SUITE(range_suite)

BOOST_AUTO_TEST_CASE(constructor_max_only_test) {
    Range r(7);
    BOOST_TEST(r.min == 0);
    BOOST_TEST(r.max == 7);
}

BOOST_AUTO_TEST_CASE(constructor_min_and_max_test) {
    Range r(7, 13);
    BOOST_TEST(r.min == 7);
    BOOST_TEST(r.max == 13);
}

BOOST_AUTO_TEST_CASE(begin_test) {
    Range r(7, 13);
    BOOST_TEST(r.begin().index == 7);
}

BOOST_AUTO_TEST_CASE(end_test) {
    Range r(7, 13);
    BOOST_TEST(r.end().index == 13);
}

BOOST_AUTO_TEST_CASE(itr_star_op_test) {
    Range r(7, 13);
    BOOST_TEST(*r.begin() == 7);
}

BOOST_AUTO_TEST_CASE(itr_plus_plus_op_test) {
    Range r(7, 13);
    RangeIterator itr = r.begin();
    BOOST_TEST(*itr == 7);
    ++itr;
    BOOST_TEST(*itr == 8);
    ++itr;
    BOOST_TEST(*itr == 9);
    ++itr;
    BOOST_TEST(*itr == 10);
    ++itr;
    BOOST_TEST(*itr == 11);
    ++itr;
    BOOST_TEST(*itr == 12);
    ++itr;
    BOOST_TEST(*itr == 13);
}

BOOST_AUTO_TEST_CASE(itr_not_equal_op_test) {
    Range r(7, 13);
    BOOST_TEST((r.begin() != r.end()));
    BOOST_TEST(!(r.begin() != r.begin()));
    RangeIterator itr = r.begin();
    ++itr;
    BOOST_TEST((r.begin() != itr));
}

BOOST_AUTO_TEST_SUITE_END()  // range_suite
BOOST_AUTO_TEST_SUITE_END()  // util_suite

}  // namespace fluoroseq
