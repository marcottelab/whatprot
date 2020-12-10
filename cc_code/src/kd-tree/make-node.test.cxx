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
#include "make-node.h"

// Standard C++ library headers:
#include <typeinfo>
#include <vector>

// Local project headers:
#include "kd-tree/internal-node.h"
#include "kd-tree/leaf-node.h"
#include "kd-tree/node.h"

namespace fluoroseq {
namespace kd_tree {

namespace {
using boost::unit_test::tolerance;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(kd_tree_suite);
BOOST_AUTO_TEST_SUITE(internal_node_suite);

BOOST_AUTO_TEST_CASE(leaf_node_test, *tolerance(TOL)) {
    int k = 2;
    int d = 3;
    vector<vector<double>> vecs(2, vector<double>(3, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 10.0;
    vecs[0][2] = 0.2;
    vecs[1][0] = 1.0;
    vecs[1][1] = 20.0;
    vecs[1][2] = 1.2;
    Node<vector<double>, vector<double>>* node =
            make_node<vector<double>, vector<double>>(k, d, &vecs[0], &vecs[2]);
    BOOST_REQUIRE(typeid(*node)
                  == typeid(LeafNode<vector<double>, vector<double>>));
    LeafNode<vector<double>, vector<double>>* leaf_node;
    leaf_node = dynamic_cast<LeafNode<vector<double>, vector<double>>*>(node);
    BOOST_TEST(leaf_node->d == d);
    BOOST_TEST(leaf_node->begin == &vecs[0]);
    BOOST_TEST(leaf_node->end == &vecs[2]);
    delete node;
}

BOOST_AUTO_TEST_CASE(internal_node_test, *tolerance(TOL)) {
    int k = 2;
    int d = 3;
    vector<vector<double>> vecs(4, vector<double>(3, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 10.0;
    vecs[0][2] = 0.2;
    vecs[1][0] = 1.0;
    vecs[1][1] = 20.0;
    vecs[1][2] = 1.2;
    vecs[2][0] = 2.0;
    vecs[2][1] = 30.0;
    vecs[2][2] = 2.2;
    vecs[3][0] = 3.0;
    vecs[3][1] = 40.0;
    vecs[3][2] = 3.2;
    Node<vector<double>, vector<double>>* node =
            make_node<vector<double>, vector<double>>(k, d, &vecs[0], &vecs[4]);
    BOOST_REQUIRE(typeid(*node)
                  == typeid(InternalNode<vector<double>, vector<double>>));
    InternalNode<vector<double>, vector<double>>* internal_node;
    internal_node =
            dynamic_cast<InternalNode<vector<double>, vector<double>>*>(node);
    BOOST_TEST(internal_node->max_left == 20.0);
    BOOST_TEST(internal_node->min_right == 30.0);
    BOOST_TEST(internal_node->split_value == 25.0);
    BOOST_TEST(internal_node->s == 1);
    BOOST_REQUIRE(typeid(*internal_node->left_child)
                  == typeid(LeafNode<vector<double>, vector<double>>));
    LeafNode<vector<double>, vector<double>>* left_child;
    left_child = dynamic_cast<LeafNode<vector<double>, vector<double>>*>(
            internal_node->left_child);
    BOOST_TEST(left_child->d == d);
    BOOST_TEST(left_child->begin == &vecs[0]);
    BOOST_TEST(left_child->end == &vecs[2]);
    vector<double> v0(3, 0);
    v0[0] = 0.0;
    v0[1] = 10.0;
    v0[2] = 0.2;
    BOOST_TEST(((left_child->begin[0] == v0) || (left_child->begin[1] == v0)));
    vector<double> v1(3, 0);
    v1[0] = 1.0;
    v1[1] = 20.0;
    v1[2] = 1.2;
    BOOST_TEST(((left_child->begin[0] == v1) || (left_child->begin[1] == v1)));
    BOOST_REQUIRE(typeid(*internal_node->right_child)
                  == typeid(LeafNode<vector<double>, vector<double>>));
    LeafNode<vector<double>, vector<double>>* right_child;
    right_child = dynamic_cast<LeafNode<vector<double>, vector<double>>*>(
            internal_node->right_child);
    BOOST_TEST(right_child->d == d);
    BOOST_TEST(right_child->begin == &vecs[2]);
    BOOST_TEST(right_child->end == &vecs[4]);
    vector<double> v2(3, 0);
    v2[0] = 2.0;
    v2[1] = 30.0;
    v2[2] = 2.2;
    BOOST_TEST(
            ((right_child->begin[0] == v2) || (right_child->begin[1] == v2)));
    vector<double> v3(3, 0);
    v3[0] = 3.0;
    v3[1] = 40.0;
    v3[2] = 3.2;
    BOOST_TEST(
            ((right_child->begin[0] == v3) || (right_child->begin[1] == v3)));
    delete node;
}

BOOST_AUTO_TEST_CASE(internal_node_reordering_test, *tolerance(TOL)) {
    int k = 2;
    int d = 3;
    vector<vector<double>> vecs(4, vector<double>(3, 0));
    vecs[0][0] = 0.0;
    vecs[0][1] = 41.0;
    vecs[0][2] = 0.2;
    vecs[1][0] = 1.0;
    vecs[1][1] = 21.0;
    vecs[1][2] = 1.2;
    vecs[2][0] = 2.0;
    vecs[2][1] = 31.0;
    vecs[2][2] = 2.2;
    vecs[3][0] = 3.0;
    vecs[3][1] = 11.0;
    vecs[3][2] = 3.2;
    Node<vector<double>, vector<double>>* node =
            make_node<vector<double>, vector<double>>(k, d, &vecs[0], &vecs[4]);
    BOOST_REQUIRE(typeid(*node)
                  == typeid(InternalNode<vector<double>, vector<double>>));
    InternalNode<vector<double>, vector<double>>* internal_node;
    internal_node =
            dynamic_cast<InternalNode<vector<double>, vector<double>>*>(node);
    BOOST_TEST(internal_node->max_left == 21.0);
    BOOST_TEST(internal_node->min_right == 31.0);
    BOOST_TEST(internal_node->split_value == 26.0);
    BOOST_TEST(internal_node->s == 1);
    BOOST_REQUIRE(typeid(*internal_node->left_child)
                  == typeid(LeafNode<vector<double>, vector<double>>));
    LeafNode<vector<double>, vector<double>>* left_child;
    left_child = dynamic_cast<LeafNode<vector<double>, vector<double>>*>(
            internal_node->left_child);
    BOOST_TEST(left_child->d == d);
    BOOST_TEST(left_child->begin == &vecs[0]);
    BOOST_TEST(left_child->end == &vecs[2]);
    vector<double> v0(3, 0);
    v0[0] = 3.0;
    v0[1] = 11.0;
    v0[2] = 3.2;
    BOOST_TEST(((left_child->begin[0] == v0) || (left_child->begin[1] == v0)));
    vector<double> v1(3, 0);
    v1[0] = 1.0;
    v1[1] = 21.0;
    v1[2] = 1.2;
    BOOST_TEST(((left_child->begin[0] == v1) || (left_child->begin[1] == v1)));
    BOOST_REQUIRE(typeid(*internal_node->right_child)
                  == typeid(LeafNode<vector<double>, vector<double>>));
    LeafNode<vector<double>, vector<double>>* right_child;
    right_child = dynamic_cast<LeafNode<vector<double>, vector<double>>*>(
            internal_node->right_child);
    BOOST_TEST(right_child->d == d);
    BOOST_TEST(right_child->begin == &vecs[2]);
    BOOST_TEST(right_child->end == &vecs[4]);
    vector<double> v2(3, 0);
    v2[0] = 2.0;
    v2[1] = 31.0;
    v2[2] = 2.2;
    BOOST_TEST(
            ((right_child->begin[0] == v2) || (right_child->begin[1] == v2)));
    vector<double> v3(3, 0);
    v3[0] = 0.0;
    v3[1] = 41.0;
    v3[2] = 0.2;
    BOOST_TEST(
            ((right_child->begin[0] == v3) || (right_child->begin[1] == v3)));
    delete node;
}

BOOST_AUTO_TEST_SUITE_END();  // internal_node_suite
BOOST_AUTO_TEST_SUITE_END();  // kd_tree_suite

}  // namespace kd_tree
}  // namespace fluoroseq
