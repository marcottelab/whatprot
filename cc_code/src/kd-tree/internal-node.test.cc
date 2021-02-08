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
#include "internal-node.h"

// Standard C++ library headers:
#include <vector>

// External headers:
#include "fakeit.hpp"

// Local project headers:
#include "kd-tree/node.h"

namespace whatprot {
namespace kd_tree {

namespace {
using boost::unit_test::tolerance;
using fakeit::_;
using fakeit::Fake;
using fakeit::Mock;
using fakeit::Verify;
using fakeit::VerifyNoOtherInvocations;
using fakeit::When;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

// Class we can use as a template parameter for E.
class Vec {
public:
    Vec(vector<double> v) : v(v), hits(1) {}
    Vec(vector<double> v, int hits) : v(v), hits(hits) {}
    double& operator[](int d) {
        return v[d];
    }
    double operator[](int d) const {
        return v[d];
    }
    vector<double> v;
    int hits;
};

BOOST_AUTO_TEST_SUITE(kd_tree_suite)
BOOST_AUTO_TEST_SUITE(internal_node_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    Mock<Node<Vec, vector<double>>> left_mock;
    Fake(Method(left_mock, search));
    Mock<Node<Vec, vector<double>>> right_mock;
    Fake(Method(right_mock, search));
    double max_left = 1.1;
    double min_right = 3.3;
    int s = 1;
    InternalNode<Vec, vector<double>> node(
            &left_mock.get(), &right_mock.get(), max_left, min_right, s);
    BOOST_TEST(node.left_child == &left_mock.get());
    BOOST_TEST(node.right_child == &right_mock.get());
    BOOST_TEST(node.max_left == 1.1);
    BOOST_TEST(node.min_right == 3.3);
    BOOST_TEST(node.split_value == (1.1 + 3.3) * 0.5);
    BOOST_TEST(node.s == 1);
}

BOOST_AUTO_TEST_CASE(search_left_only_test, *tolerance(TOL)) {
    Mock<Node<Vec, vector<double>>> left_mock;
    When(Method(left_mock, search))
            .Do([](const vector<double>& query, KBest<Vec>* k_best) {
                k_best->kth_dist_sq = 1.0 * 1.0;
            });
    Mock<Node<Vec, vector<double>>> right_mock;
    Fake(Method(right_mock, search));
    double max_left = 1.1;
    double min_right = 3.3;
    int s = 1;
    InternalNode<Vec, vector<double>> node(
            &left_mock.get(), &right_mock.get(), max_left, min_right, s);
    vector<double> query(3, 0);
    query[0] = 30;
    query[1] = 2.0;
    query[2] = 40;
    KBest<Vec> k_best(0);
    node.search(query, &k_best);
    Verify(Method(left_mock, search).Using(query, _)).Exactly(1);
    VerifyNoOtherInvocations(left_mock);
    VerifyNoOtherInvocations(right_mock);
}

BOOST_AUTO_TEST_CASE(search_left_check_using_squared_distances_test,
                     *tolerance(TOL)) {
    Mock<Node<Vec, vector<double>>> left_mock;
    When(Method(left_mock, search))
            .Do([](const vector<double>& query, KBest<Vec>* k_best) {
                // This value is carefully chosen so that it is just small
                // enough to not trigger the second search, but big enough so
                // that if the calculation is wrong, not accounting for the fact
                // that this is a squared distance, it will trigger the second
                // search.
                k_best->kth_dist_sq = 1.2 * 1.2;
            });
    Mock<Node<Vec, vector<double>>> right_mock;
    Fake(Method(right_mock, search));
    double max_left = 1.1;
    double min_right = 3.3;
    int s = 1;
    InternalNode<Vec, vector<double>> node(
            &left_mock.get(), &right_mock.get(), max_left, min_right, s);
    vector<double> query(3, 0);
    query[0] = 30;
    query[1] = 2.0;
    query[2] = 40;
    KBest<Vec> k_best(0);
    node.search(query, &k_best);
    Verify(Method(left_mock, search).Using(query, _)).Exactly(1);
    VerifyNoOtherInvocations(left_mock);
    VerifyNoOtherInvocations(right_mock);
}

BOOST_AUTO_TEST_CASE(search_left_then_right_test, *tolerance(TOL)) {
    Mock<Node<Vec, vector<double>>> left_mock;
    When(Method(left_mock, search))
            .Do([](const vector<double>& query, KBest<Vec>* k_best) {
                k_best->kth_dist_sq = 2.0 * 2.0;
            });
    Mock<Node<Vec, vector<double>>> right_mock;
    Fake(Method(right_mock, search));
    double max_left = 1.1;
    double min_right = 3.3;
    int s = 1;
    InternalNode<Vec, vector<double>> node(
            &left_mock.get(), &right_mock.get(), max_left, min_right, s);
    vector<double> query(3, 0);
    query[0] = 30;
    query[1] = 2.0;
    query[2] = 40;
    KBest<Vec> k_best(0);
    node.search(query, &k_best);
    Verify(Method(left_mock, search).Using(query, _)
           + Method(right_mock, search).Using(query, _))
            .Exactly(1);
    VerifyNoOtherInvocations(left_mock);
    VerifyNoOtherInvocations(right_mock);
}

BOOST_AUTO_TEST_CASE(search_right_only_test, *tolerance(TOL)) {
    Mock<Node<Vec, vector<double>>> left_mock;
    Fake(Method(left_mock, search));
    Mock<Node<Vec, vector<double>>> right_mock;
    When(Method(right_mock, search))
            .Do([](const vector<double>& query, KBest<Vec>* k_best) {
                k_best->kth_dist_sq = 1.0 * 1.0;
            });
    double max_left = 1.1;
    double min_right = 3.3;
    int s = 1;
    InternalNode<Vec, vector<double>> node(
            &left_mock.get(), &right_mock.get(), max_left, min_right, s);
    vector<double> query(3, 0);
    query[0] = 30;
    query[1] = 3.0;
    query[2] = 40;
    KBest<Vec> k_best(0);
    node.search(query, &k_best);
    Verify(Method(right_mock, search).Using(query, _)).Exactly(1);
    VerifyNoOtherInvocations(right_mock);
    VerifyNoOtherInvocations(left_mock);
}

BOOST_AUTO_TEST_CASE(search_right_check_using_squared_distances_test,
                     *tolerance(TOL)) {
    Mock<Node<Vec, vector<double>>> left_mock;
    Fake(Method(left_mock, search));
    Mock<Node<Vec, vector<double>>> right_mock;
    When(Method(right_mock, search))
            .Do([](const vector<double>& query, KBest<Vec>* k_best) {
                // This value is carefully chosen so that it is just small
                // enough to not trigger the second search, but big enough so
                // that if the calculation is wrong, not accounting for the fact
                // that this is a squared distance, it will trigger the second
                // search.
                k_best->kth_dist_sq = 1.2 * 1.2;
            });
    double max_left = 1.1;
    double min_right = 3.3;
    int s = 1;
    InternalNode<Vec, vector<double>> node(
            &left_mock.get(), &right_mock.get(), max_left, min_right, s);
    vector<double> query(3, 0);
    query[0] = 30;
    query[1] = 3.0;
    query[2] = 40;
    KBest<Vec> k_best(0);
    node.search(query, &k_best);
    Verify(Method(right_mock, search).Using(query, _)).Exactly(1);
    VerifyNoOtherInvocations(right_mock);
    VerifyNoOtherInvocations(left_mock);
}

BOOST_AUTO_TEST_CASE(search_right_then_left_test, *tolerance(TOL)) {
    Mock<Node<Vec, vector<double>>> left_mock;
    Fake(Method(left_mock, search));
    Mock<Node<Vec, vector<double>>> right_mock;
    When(Method(right_mock, search))
            .Do([](const vector<double>& query, KBest<Vec>* k_best) {
                k_best->kth_dist_sq = 2.0 * 2.0;
            });
    double max_left = 1.1;
    double min_right = 3.3;
    int s = 1;
    InternalNode<Vec, vector<double>> node(
            &left_mock.get(), &right_mock.get(), max_left, min_right, s);
    vector<double> query(3, 0);
    query[0] = 30;
    query[1] = 3.0;
    query[2] = 40;
    KBest<Vec> k_best(0);
    node.search(query, &k_best);
    Verify(Method(right_mock, search).Using(query, _)
           + Method(left_mock, search).Using(query, _))
            .Exactly(1);
    VerifyNoOtherInvocations(right_mock);
    VerifyNoOtherInvocations(left_mock);
}

BOOST_AUTO_TEST_SUITE_END()  // internal_node_suite
BOOST_AUTO_TEST_SUITE_END()  // kd_tree_suite

}  // namespace kd_tree
}  // namespace whatprot
