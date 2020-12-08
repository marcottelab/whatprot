/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin *
* Department: Oden Institute and Institute for Cellular and Molecular Biology
*
* PI: Edward Marcotte *
* Project: Protein Fluorosequencing *
\******************************************************************************/

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// File under test:
#include "kd-tree.h"

// Standard C++ library headers:
#include <algorithm>
#include <vector>

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
using std::move;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(kd_tree_suite);
BOOST_AUTO_TEST_SUITE(kd_tree_suite);

BOOST_AUTO_TEST_CASE(big_test, *tolerance(TOL)) {
    int k = 4;
    int d = 2;
    // We'll put in a 4x3 grid of points, all on integers.
    vector<vector<double>> vecs(12, vector<double>(2, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 0.0;
    vecs[1][0] = 0.0;
    vecs[1][1] = 1.0;
    vecs[2][0] = 0.0;
    vecs[2][1] = 2.0;
    vecs[3][0] = 0.0;
    vecs[3][1] = 3.0;
    vecs[4][0] = 1.0;
    vecs[4][1] = 0.0;
    vecs[5][0] = 1.0;
    vecs[5][1] = 1.0;
    vecs[6][0] = 1.0;
    vecs[6][1] = 2.0;
    vecs[7][0] = 1.0;
    vecs[7][1] = 3.0;
    vecs[8][0] = 2.0;
    vecs[8][1] = 0.0;
    vecs[9][0] = 2.0;
    vecs[9][1] = 1.0;
    vecs[10][0] = 2.0;
    vecs[10][1] = 2.0;
    vecs[11][0] = 2.0;
    vecs[11][1] = 3.0;
    KDTree<vector<double>, vector<double>> kdt(k, d, move(vecs));
    vector<double> query(2, 0);
    query[0] = .9;
    query[1] = .8;
    vector<vector<double>*> k_nearest;
    vector<double> dists_sq;
    kdt.search(query, &k_nearest, &dists_sq);
    vector<double> n4(2, 0);
    n4[0] = 2.0;
    n4[1] = 1.0;
    BOOST_TEST(*(k_nearest[0]) == n4);
    BOOST_TEST(dists_sq[0]
               == (2.0 - 0.9) * (2.0 - 0.9) + (1.0 - 0.8) * (1.0 - 0.8));
    vector<double> n3(2, 0);
    n3[0] = 0.0;
    n3[1] = 1.0;
    BOOST_TEST(*(k_nearest[1]) == n3);
    BOOST_TEST(dists_sq[1]
               == (0.0 - 0.9) * (0.0 - 0.9) + (1.0 - 0.8) * (1.0 - 0.8));
    vector<double> n2(2, 0);
    n2[0] = 1.0;
    n2[1] = 0.0;
    BOOST_TEST(*(k_nearest[2]) == n2);
    BOOST_TEST(dists_sq[2]
               == (1.0 - 0.9) * (1.0 - 0.9) + (0.0 - 0.8) * (0.0 - 0.8));
    vector<double> n1(2, 0);
    n1[0] = 1.0;
    n1[1] = 1.0;
    BOOST_TEST(*(k_nearest[3]) == n1);
    BOOST_TEST(dists_sq[3]
               == (1.0 - 0.9) * (1.0 - 0.9) + (1.0 - 0.8) * (1.0 - 0.8));
}

BOOST_AUTO_TEST_SUITE_END();  // kd_tree_suite
BOOST_AUTO_TEST_SUITE_END();  // kd_tree_suite

}  // namespace fluoroseq
