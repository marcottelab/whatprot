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
#include "k-best.h"

// Standard C++ library headers:
#include <cfloat>
#include <string>
#include <vector>

namespace fluoroseq {
namespace kd_tree {

namespace {
using boost::unit_test::tolerance;
using std::string;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(kd_tree_suite);
BOOST_AUTO_TEST_SUITE(k_best_suite);

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    KBest<string> kb(3);
    BOOST_TEST(kb.k == 3);
    BOOST_TEST(kb.kth_distance == DBL_MAX);
}

BOOST_AUTO_TEST_CASE(empty_test, *tolerance(TOL)) {
    KBest<string> kb(3);
    vector<string*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_TEST(v.size() == 0);
    BOOST_TEST(dists_sq.size() == 0);
}

BOOST_AUTO_TEST_CASE(under_capacity_test, *tolerance(TOL)) {
    KBest<string> kb(3);
    string s_one = "one point zero";
    kb.consider(1.0, &s_one);
    BOOST_TEST(kb.kth_distance == DBL_MAX);
    vector<string*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_TEST(v.size() == 1);
    BOOST_TEST(*v[0] == "one point zero");
    BOOST_TEST(dists_sq[0] == 1.0);
}

BOOST_AUTO_TEST_CASE(at_capacity_test, *tolerance(TOL)) {
    KBest<string> kb(3);
    string s_one = "one point zero";
    kb.consider(1.0, &s_one);
    string s_three = "three point zero";
    kb.consider(3.0, &s_three);
    string s_two = "two point zero";
    kb.consider(2.0, &s_two);
    BOOST_TEST(kb.kth_distance == 3.0);
    vector<string*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_TEST(v.size() == 3);
    BOOST_TEST(dists_sq.size() == 3);
    BOOST_TEST(*v[0] == "three point zero");
    BOOST_TEST(dists_sq[0] == 3.0);
    BOOST_TEST(*v[1] == "two point zero");
    BOOST_TEST(dists_sq[1] == 2.0);
    BOOST_TEST(*v[2] == "one point zero");
    BOOST_TEST(dists_sq[2] == 1.0);
}

BOOST_AUTO_TEST_CASE(over_capacity_test, *tolerance(TOL)) {
    KBest<string> kb(3);
    string s_one = "one point zero";
    kb.consider(1.0, &s_one);
    string s_three = "three point zero";
    kb.consider(3.0, &s_three);
    string s_two = "two point zero";
    kb.consider(2.0, &s_two);
    string s_first_extra = "four point four";
    kb.consider(4.4, &s_first_extra);
    string s_second_extra = "two point two";
    kb.consider(2.2, &s_second_extra);
    BOOST_TEST(kb.kth_distance == 2.2);
    vector<string*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_TEST(v.size() == 3);
    BOOST_TEST(dists_sq.size() == 3);
    BOOST_TEST(*v[0] == "two point two");
    BOOST_TEST(dists_sq[0] == 2.2);
    BOOST_TEST(*v[1] == "two point zero");
    BOOST_TEST(dists_sq[1] == 2.0);
    BOOST_TEST(*v[2] == "one point zero");
    BOOST_TEST(dists_sq[2] == 1.0);
}

BOOST_AUTO_TEST_SUITE_END();  // k_best_suite
BOOST_AUTO_TEST_SUITE_END();  // kd_tree_suite

}  // namespace kd_tree
}  // namespace fluoroseq
