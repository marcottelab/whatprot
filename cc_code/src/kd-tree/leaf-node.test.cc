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
    // Needed to check arguments with FakeIt.
    bool operator==(const Vec& other) const {
        return (v == other.v) && (hits == other.hits);
    }
    vector<double> v;
    int hits;
};

BOOST_AUTO_TEST_SUITE(kd_tree_suite)
BOOST_AUTO_TEST_SUITE(leaf_node_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    int d = 3;
    vector<Vec> vecs(2, Vec(vector<double>(3, 0)));
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    BOOST_TEST(leaf.d == d);
    BOOST_TEST(leaf.begin == &vecs[0]);
    BOOST_TEST(leaf.end == &vecs[2]);
}

BOOST_AUTO_TEST_CASE(consider_success_test, *tolerance(TOL)) {
    int d = 3;
    vector<Vec> vecs(2, Vec(vector<double>(3, 0)));
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    vector<double> query(3, 0);
    query[0] = 1.0;
    query[1] = 1.1;
    query[2] = 1.2;
    Vec entry(vector<double>(3, 0));
    entry[0] = 2.00;
    entry[1] = 2.11;
    entry[2] = 2.22;
    Mock<KBest<Vec>> k_best_mock;
    Fake(Method(k_best_mock, insert));
    KBest<Vec>* k_best = &k_best_mock.get();
    k_best->kth_dist_sq = 1e9;
    leaf.consider(query, &entry, k_best);
    double dist_sq = 1.0 * 1.0 + 1.01 * 1.01 + 1.02 * 1.02;
    Verify(Method(k_best_mock, insert).Using(Close(dist_sq, TOL), Ptr(entry)));
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_CASE(consider_failure_test, *tolerance(TOL)) {
    int d = 3;
    vector<Vec> vecs(2, Vec(vector<double>(3, 0)));
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    vector<double> query(3, 0);
    query[0] = 1.0;
    query[1] = 1.1;
    query[2] = 1.2;
    Vec entry(vector<double>(3, 0));
    entry[0] = 2.00;
    entry[1] = 2.11;
    entry[2] = 2.22;
    Mock<KBest<Vec>> k_best_mock;
    Fake(Method(k_best_mock, insert));
    KBest<Vec>* k_best = &k_best_mock.get();
    k_best->kth_dist_sq = 1e-2;
    leaf.consider(query, &entry, k_best);
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_CASE(consider_success_big_d_test, *tolerance(TOL)) {
    int d = 6;
    vector<Vec> vecs(2, Vec(vector<double>(6, 0)));
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    vector<double> query(6, 0);
    query[0] = 1.0;
    query[1] = 1.1;
    query[2] = 1.2;
    query[3] = 1.3;
    query[4] = 1.4;
    query[5] = 1.5;
    Vec entry(vector<double>(6, 0));
    entry[0] = 2.00;
    entry[1] = 2.11;
    entry[2] = 2.22;
    entry[3] = 2.33;
    entry[4] = 2.44;
    entry[5] = 2.55;
    Mock<KBest<Vec>> k_best_mock;
    Fake(Method(k_best_mock, insert));
    KBest<Vec>* k_best = &k_best_mock.get();
    k_best->kth_dist_sq = 1e9;
    leaf.consider(query, &entry, k_best);
    double dist_sq = 1.0 * 1.0 + 1.01 * 1.01 + 1.02 * 1.02 + 1.03 * 1.03
                     + 1.04 * 1.04 + 1.05 * 1.05;
    Verify(Method(k_best_mock, insert).Using(Close(dist_sq, TOL), Ptr(entry)));
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_CASE(consider_failure_big_d_test, *tolerance(TOL)) {
    int d = 6;
    vector<Vec> vecs(2, Vec(vector<double>(6, 0)));
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    vector<double> query(6, 0);
    query[0] = 1.0;
    query[1] = 1.1;
    query[2] = 1.2;
    query[3] = 1.3;
    query[4] = 1.4;
    query[5] = 1.5;
    Vec entry(vector<double>(6, 0));
    entry[0] = 2.00;
    entry[1] = 2.11;
    entry[2] = 2.22;
    entry[3] = 2.33;
    entry[4] = 2.44;
    entry[5] = 2.55;
    Mock<KBest<Vec>> k_best_mock;
    Fake(Method(k_best_mock, insert));
    KBest<Vec>* k_best = &k_best_mock.get();
    k_best->kth_dist_sq = 1e-2;
    leaf.consider(query, &entry, k_best);
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_CASE(consider_barely_failure_big_d_test, *tolerance(TOL)) {
    int d = 6;
    vector<Vec> vecs(2, Vec(vector<double>(6, 0)));
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    vector<double> query(6, 0);
    query[0] = 1.0;
    query[1] = 1.1;
    query[2] = 1.2;
    query[3] = 1.3;
    query[4] = 1.4;
    query[5] = 1.5;
    Vec entry(vector<double>(6, 0));
    entry[0] = 2.00;
    entry[1] = 2.11;
    entry[2] = 2.22;
    entry[3] = 2.33;
    entry[4] = 2.44;
    entry[5] = 2.55;
    Mock<KBest<Vec>> k_best_mock;
    Fake(Method(k_best_mock, insert));
    KBest<Vec>* k_best = &k_best_mock.get();
    k_best->kth_dist_sq =
            1.0 * 1.0 + 1.01 * 1.01 + 1.02 * 1.02 + 1.03 * 1.03 - 1e-7;
    leaf.consider(query, &entry, k_best);
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_CASE(search_success_test, *tolerance(TOL)) {
    int d = 3;
    vector<Vec> vecs(2, Vec(vector<double>(3, 0)));
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
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    Mock<KBest<Vec>> k_best_mock;
    Fake(Method(k_best_mock, insert));
    KBest<Vec>* k_best = &k_best_mock.get();
    k_best->kth_dist_sq = 1e9;
    leaf.search(query, k_best);
    double dst1 = 0.0;
    dst1 += (0.2 - 0.0) * (0.2 - 0.0);
    dst1 += (0.3 - 0.1) * (0.3 - 0.1);
    dst1 += (0.5 - 0.2) * (0.5 - 0.2);
    Vec v1(vector<double>(3, 0));
    v1[0] = 0.0;
    v1[1] = 0.1;
    v1[2] = 0.2;
    Verify(Method(k_best_mock, insert).Using(Close(dst1, TOL), Ptr(v1)))
            .Exactly(1);
    double dst2 = 0.0;
    dst2 += (0.2 - 1.0) * (0.2 - 1.0);
    dst2 += (0.3 - 1.1) * (0.3 - 1.1);
    dst2 += (0.5 - 1.2) * (0.5 - 1.2);
    Vec v2(vector<double>(3, 0));
    v2[0] = 1.0;
    v2[1] = 1.1;
    v2[2] = 1.2;
    Verify(Method(k_best_mock, insert).Using(Close(dst2, TOL), Ptr(v2)))
            .Exactly(1);
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_CASE(search_failure_test, *tolerance(TOL)) {
    int d = 3;
    vector<Vec> vecs(2, Vec(vector<double>(3, 0)));
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
    LeafNode<Vec, vector<double>> leaf(d, &vecs[0], &vecs[2]);
    Mock<KBest<Vec>> k_best_mock;
    Fake(Method(k_best_mock, insert));
    KBest<Vec>* k_best = &k_best_mock.get();
    k_best->kth_dist_sq = 1e-2;
    leaf.search(query, k_best);
    VerifyNoOtherInvocations(k_best_mock);
}

BOOST_AUTO_TEST_SUITE_END()  // leaf_node_suite
BOOST_AUTO_TEST_SUITE_END()  // kd_tree_suite

}  // namespace kd_tree
}  // namespace fluoroseq
