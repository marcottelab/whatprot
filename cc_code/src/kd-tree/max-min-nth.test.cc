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
#include "max-min-nth.h"

// Standard C++ library headers:
#include <vector>

namespace whatprot {
namespace kd_tree {

namespace {
using boost::unit_test::tolerance;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(kd_tree_suite)
BOOST_AUTO_TEST_SUITE(max_min_nth_suite)

BOOST_AUTO_TEST_CASE(partition_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 99.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 11.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 50.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 88.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 22.0;
    int num_eq_to_pivot;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth =
            partition(&vecs[0], &vecs[5], 1, nth, &num_eq_to_pivot);
    BOOST_TEST(num_eq_to_pivot == 0);
    BOOST_TEST(new_nth == &vecs[2]);
    vector<double> low_a(2, 0);
    low_a[0] = 1.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 4.0;
    low_b[1] = 22.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> pivot(2, 0);
    pivot[0] = 2.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[2] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 3.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 0.0;
    high_b[1] = 99.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(partition_with_duplicates_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 50.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 99.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 50.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 22.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 50.0;
    int num_eq_to_pivot;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth =
            partition(&vecs[0], &vecs[5], 1, nth, &num_eq_to_pivot);
    BOOST_TEST(num_eq_to_pivot == 2);
    BOOST_TEST(new_nth == &vecs[3]);
    vector<double> low_a(2, 0);
    low_a[0] = 0.0;
    low_a[1] = 50.0;
    BOOST_TEST(
            ((vecs[0] == low_a) || (vecs[1] == low_a) || (vecs[2] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 3.0;
    low_b[1] = 22.0;
    BOOST_TEST(
            ((vecs[0] == low_b) || (vecs[1] == low_b) || (vecs[2] == low_b)));
    vector<double> low_c(2, 0);
    low_c[0] = 4.0;
    low_c[1] = 50.0;
    BOOST_TEST(
            ((vecs[0] == low_c) || (vecs[1] == low_c) || (vecs[2] == low_c)));
    vector<double> pivot(2, 0);
    pivot[0] = 2.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[3] == pivot);
    vector<double> high(2, 0);
    high[0] = 1.0;
    high[1] = 99.0;
    BOOST_TEST(vecs[4] == high);
}

BOOST_AUTO_TEST_CASE(partition_alternate_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 99.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 11.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 22.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 88.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 50.0;
    partition_alternate(&vecs[0], &vecs[5], 1);
    vector<double> low_a(2, 0);
    low_a[0] = 1.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 2.0;
    low_b[1] = 22.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> pivot(2, 0);
    pivot[0] = 4.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[2] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 3.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 0.0;
    high_b[1] = 99.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(partition_alternate_with_duplicates_test,
                     *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 50.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 99.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 50.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 22.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 50.0;
    partition_alternate(&vecs[0], &vecs[5], 1);
    vector<double> low(2, 0);
    low[0] = 3.0;
    low[1] = 22.0;
    BOOST_TEST(vecs[0] == low);
    vector<double> pivot(2, 0);
    pivot[0] = 4.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[1] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 0.0;
    high_a[1] = 50.0;
    BOOST_TEST(((vecs[2] == high_a) || (vecs[3] == high_a)
                || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 2.0;
    high_b[1] = 50.0;
    BOOST_TEST(((vecs[2] == high_b) || (vecs[3] == high_b)
                || (vecs[4] == high_b)));
    vector<double> high_c(2, 0);
    high_c[0] = 1.0;
    high_c[1] = 99.0;
    BOOST_TEST(((vecs[2] == high_c) || (vecs[3] == high_c)
                || (vecs[4] == high_c)));
}

BOOST_AUTO_TEST_CASE(nth_element_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 99.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 11.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 50.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 88.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 22.0;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth = nth_element(&vecs[0], &vecs[5], 1, nth);
    BOOST_TEST(new_nth == &vecs[2]);
    vector<double> low_a(2, 0);
    low_a[0] = 1.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 4.0;
    low_b[1] = 22.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> pivot(2, 0);
    pivot[0] = 2.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[2] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 3.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 0.0;
    high_b[1] = 99.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(nth_element_recursive_on_left_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 99.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 11.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 88.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 50.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 22.0;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth = nth_element(&vecs[0], &vecs[5], 1, nth);
    BOOST_TEST(new_nth == &vecs[2]);
    vector<double> low_a(2, 0);
    low_a[0] = 1.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 4.0;
    low_b[1] = 22.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> pivot(2, 0);
    pivot[0] = 3.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[2] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 2.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 0.0;
    high_b[1] = 99.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(nth_element_recursive_on_left_with_duplicate_test,
                     *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 88.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 11.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 88.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 50.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 22.0;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth = nth_element(&vecs[0], &vecs[5], 1, nth);
    BOOST_TEST(new_nth == &vecs[2]);
    vector<double> low_a(2, 0);
    low_a[0] = 1.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 4.0;
    low_b[1] = 22.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> pivot(2, 0);
    pivot[0] = 3.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[2] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 2.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 0.0;
    high_b[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(nth_element_recursive_on_right_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 99.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 50.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 11.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 88.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 22.0;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth = nth_element(&vecs[0], &vecs[5], 1, nth);
    BOOST_TEST(new_nth == &vecs[2]);
    vector<double> low_a(2, 0);
    low_a[0] = 2.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 4.0;
    low_b[1] = 22.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> pivot(2, 0);
    pivot[0] = 1.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[2] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 3.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 0.0;
    high_b[1] = 99.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(nth_element_recursive_on_right_with_duplicate_test,
                     *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 99.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 50.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 11.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 88.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 11.0;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth = nth_element(&vecs[0], &vecs[5], 1, nth);
    BOOST_TEST(new_nth == &vecs[2]);
    vector<double> low_a(2, 0);
    low_a[0] = 2.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 4.0;
    low_b[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> pivot(2, 0);
    pivot[0] = 1.0;
    pivot[1] = 50.0;
    BOOST_TEST(vecs[2] == pivot);
    vector<double> high_a(2, 0);
    high_a[0] = 3.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 0.0;
    high_b[1] = 99.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(nth_element_duplicates_to_the_left_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 50.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 99.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 50.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 22.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 88.0;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth = nth_element(&vecs[0], &vecs[5], 1, nth);
    BOOST_TEST(new_nth == &vecs[3]);
    vector<double> low_a(2, 0);
    low_a[0] = 0.0;
    low_a[1] = 50.0;
    BOOST_TEST(
            ((vecs[0] == low_a) || (vecs[1] == low_a) || (vecs[2] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 3.0;
    low_b[1] = 22.0;
    BOOST_TEST(
            ((vecs[0] == low_b) || (vecs[1] == low_b) || (vecs[2] == low_b)));
    vector<double> low_c(2, 0);
    low_c[0] = 2.0;
    low_c[1] = 50.0;
    BOOST_TEST(
            ((vecs[0] == low_c) || (vecs[1] == low_c) || (vecs[2] == low_c)));
    vector<double> high_a(2, 0);
    high_a[0] = 4.0;
    high_a[1] = 88.0;
    BOOST_TEST(((vecs[3] == high_a) || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 1.0;
    high_b[1] = 99.0;
    BOOST_TEST(((vecs[3] == high_b) || (vecs[4] == high_b)));
}

BOOST_AUTO_TEST_CASE(nth_element_duplicates_to_the_right_test,
                     *tolerance(TOL)) {
    vector<vector<double>> vecs(5, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 11.0;
    vecs[1][0] = 1.0;
    vecs[1][1] = 99.0;
    vecs[2][0] = 2.0;
    vecs[2][1] = 50.0;
    vecs[3][0] = 3.0;
    vecs[3][1] = 22.0;
    vecs[4][0] = 4.0;
    vecs[4][1] = 50.0;
    vector<double>* nth = &vecs[2];
    vector<double>* new_nth = nth_element(&vecs[0], &vecs[5], 1, nth);
    BOOST_TEST(new_nth == &vecs[2]);
    vector<double> low_a(2, 0);
    low_a[0] = 0.0;
    low_a[1] = 11.0;
    BOOST_TEST(((vecs[0] == low_a) || (vecs[1] == low_a)));
    vector<double> low_b(2, 0);
    low_b[0] = 3.0;
    low_b[1] = 22.0;
    BOOST_TEST(((vecs[0] == low_b) || (vecs[1] == low_b)));
    vector<double> high_a(2, 0);
    high_a[0] = 2.0;
    high_a[1] = 50.0;
    BOOST_TEST(((vecs[2] == high_a) || (vecs[3] == high_a)
                || (vecs[4] == high_a)));
    vector<double> high_b(2, 0);
    high_b[0] = 4.0;
    high_b[1] = 50.0;
    BOOST_TEST(((vecs[2] == high_b) || (vecs[3] == high_b)
                || (vecs[4] == high_b)));
    vector<double> high_c(2, 0);
    high_c[0] = 1.0;
    high_c[1] = 99.0;
    BOOST_TEST(((vecs[2] == high_c) || (vecs[3] == high_c)
                || (vecs[4] == high_c)));
}

BOOST_AUTO_TEST_CASE(max_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(4, vector<double>(3, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 0.1;
    vecs[0][2] = 0.2;
    vecs[1][0] = 1.0;
    vecs[1][1] = 1.1;
    vecs[1][2] = 1.2;
    vecs[2][0] = 2.0;
    vecs[2][1] = 2.1;
    vecs[2][2] = 2.2;
    vecs[3][0] = 3.0;
    vecs[3][1] = -3.1;
    vecs[3][2] = 3.2;
    double max = max_element(&vecs[0], &vecs[4], 1);
    BOOST_TEST(max == 2.1);
}

BOOST_AUTO_TEST_CASE(min_test, *tolerance(TOL)) {
    vector<vector<double>> vecs(4, vector<double>(3, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 0.1;
    vecs[0][2] = 0.2;
    vecs[1][0] = 1.0;
    vecs[1][1] = -1.1;
    vecs[1][2] = 1.2;
    vecs[2][0] = 2.0;
    vecs[2][1] = 2.1;
    vecs[2][2] = 2.2;
    vecs[3][0] = 3.0;
    vecs[3][1] = 3.1;
    vecs[3][2] = 3.2;
    double min = min_element(&vecs[0], &vecs[4], 1);
    BOOST_TEST(min == -1.1);
}

BOOST_AUTO_TEST_SUITE_END()  // max_min_nth_suite
BOOST_AUTO_TEST_SUITE_END()  // kd_tree_suite

}  // namespace kd_tree
}  // namespace whatprot
