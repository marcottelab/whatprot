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
#include "leaf-node.h"

// Standard C++ library headers:
#include <vector>

// External headers:
#include "fakeit.hpp"

// Local project headers:
#include "test-util/fakeit.h"  // in test directory

namespace fluoroseq {
namespace kd_tree {

namespace {
using boost::unit_test::tolerance;
using fakeit::_;
using fakeit::Fake;
using fakeit::Mock;
using fakeit::Verify;
using fakeit::VerifyNoOtherInvocations;
using fluoroseq::test_util::Close;
using fluoroseq::test_util::Ptr;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(kd_tree_suite);
BOOST_AUTO_TEST_SUITE(leaf_node_suite);

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    int d = 3;
    vector<vector<double>> vecs(2, vector<double>(3, 0));
    LeafNode<vector<double>, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    BOOST_TEST(leaf.d == d);
    BOOST_TEST(leaf.begin == &vecs[0]);
    BOOST_TEST(leaf.end == &vecs[2]);
}

BOOST_AUTO_TEST_CASE(distance_test, *tolerance(TOL)) {
    int d = 3;
    vector<vector<double>> vecs(2, vector<double>(3, 0));
    LeafNode<vector<double>, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    vector<double> t1(3, 0);
    t1[0] = 1.0;
    t1[1] = 1.1;
    t1[2] = 1.2;
    vector<double> t2(3, 0);
    t2[0] = 2.00;
    t2[1] = 2.11;
    t2[2] = 2.22;
    BOOST_TEST(leaf.distance(t1, t2) == 1.0 * 1.0 + 1.01 * 1.01 + 1.02 * 1.02);
}

BOOST_AUTO_TEST_CASE(search_test, *tolerance(TOL)) {
    int d = 3;
    vector<vector<double>> vecs(2, vector<double>(3, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 0.1;
    vecs[0][2] = 0.2;
    vecs[1][0] = 1.0;
    vecs[1][1] = 1.1;
    vecs[1][2] = 1.2;
    vector<double> query(3, 0);
    query[0] = 0.2;
    query[1] = 0.3;
    query[2] = 0.5;
    LeafNode<vector<double>, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    Mock<KBest<vector<double>>> k_best_mock;
    Fake(Method(k_best_mock, consider));
    leaf.search(query, &k_best_mock.get());
    double dst1 = 0.0;
    dst1 += (0.2 - 0.0) * (0.2 - 0.0);
    dst1 += (0.3 - 0.1) * (0.3 - 0.1);
    dst1 += (0.5 - 0.2) * (0.5 - 0.2);
    vector<double> v1(3, 0);
    v1[0] = 0.0;
    v1[1] = 0.1;
    v1[2] = 0.2;
    Verify(Method(k_best_mock, consider).Using(Close(dst1, TOL), Ptr(v1)))
            .Exactly(1);
    double dst2 = 0.0;
    dst2 += (0.2 - 1.0) * (0.2 - 1.0);
    dst2 += (0.3 - 1.1) * (0.3 - 1.1);
    dst2 += (0.5 - 1.2) * (0.5 - 1.2);
    vector<double> v2(3, 0);
    v2[0] = 1.0;
    v2[1] = 1.1;
    v2[2] = 1.2;
    Verify(Method(k_best_mock, consider).Using(Close(dst2, TOL), Ptr(v2)))
            .Exactly(1);
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_SUITE_END();  // leaf_node_suite
BOOST_AUTO_TEST_SUITE_END();  // kd_tree_suite

}  // namespace kd_tree
}  // namespace fluoroseq
