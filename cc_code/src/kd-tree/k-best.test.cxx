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

// Class we can use as a template parameter for E. Note that we are missing the
// [] operator. This is OK. It's not needed by KBest<E>.
class Str {
public:
    Str(string s) : s(s), hits(1) {}
    Str(string s, int hits) : s(s), hits(hits) {}
    int hits;
    string s;
};

BOOST_AUTO_TEST_SUITE(kd_tree_suite);
BOOST_AUTO_TEST_SUITE(k_best_suite);

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    BOOST_TEST(kb.k == 3);
    BOOST_TEST(kb.kth_dist_sq == DBL_MAX);
    BOOST_TEST(kb.hits == 0);
}

BOOST_AUTO_TEST_CASE(empty_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    vector<Str*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_REQUIRE(v.size() == 0);
    BOOST_TEST(dists_sq.size() == 0);
}

BOOST_AUTO_TEST_CASE(under_capacity_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    Str s_one("one point zero");
    kb.insert(1.0, &s_one);
    BOOST_TEST(kb.kth_dist_sq == DBL_MAX);
    BOOST_TEST(kb.hits == 1);
    vector<Str*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_REQUIRE(v.size() == 1);
    BOOST_TEST(v[0]->s == "one point zero");
    BOOST_TEST(v[0]->hits == 1);
    BOOST_TEST(dists_sq[0] == 1.0);
}

BOOST_AUTO_TEST_CASE(under_capacity_hits_gt1_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    Str s_one("one point zero", 2);
    kb.insert(1.0, &s_one);
    BOOST_TEST(kb.kth_dist_sq == DBL_MAX);
    BOOST_TEST(kb.hits == 2);
    vector<Str*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_REQUIRE(v.size() == 1);
    BOOST_TEST(v[0]->s == "one point zero");
    BOOST_TEST(v[0]->hits == 2);
    BOOST_TEST(dists_sq[0] == 1.0);
}

BOOST_AUTO_TEST_CASE(at_capacity_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    Str s_one("one point zero");
    kb.insert(1.0, &s_one);
    Str s_three("three point zero");
    kb.insert(3.0, &s_three);
    Str s_two("two point zero");
    kb.insert(2.0, &s_two);
    BOOST_TEST(kb.kth_dist_sq == 3.0);
    BOOST_TEST(kb.hits == 3);
    vector<Str*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_REQUIRE(v.size() == 3);
    BOOST_TEST(dists_sq.size() == 3);
    BOOST_TEST(v[0]->s == "three point zero");
    BOOST_TEST(v[0]->hits == 1);
    BOOST_TEST(dists_sq[0] == 3.0);
    BOOST_TEST(v[1]->s == "two point zero");
    BOOST_TEST(v[1]->hits == 1);
    BOOST_TEST(dists_sq[1] == 2.0);
    BOOST_TEST(v[2]->s == "one point zero");
    BOOST_TEST(v[2]->hits == 1);
    BOOST_TEST(dists_sq[2] == 1.0);
}

BOOST_AUTO_TEST_CASE(at_capacity_hits_gt1_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    Str s_one("one point zero", 1);
    kb.insert(1.0, &s_one);
    Str s_two("two point zero", 2);
    kb.insert(2.0, &s_two);
    BOOST_TEST(kb.kth_dist_sq == 2.0);
    BOOST_TEST(kb.hits == 3);
    vector<Str*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_REQUIRE(v.size() == 2);
    BOOST_TEST(dists_sq.size() == 2);
    BOOST_TEST(v[0]->s == "two point zero");
    BOOST_TEST(v[0]->hits == 2);
    BOOST_TEST(dists_sq[0] == 2.0);
    BOOST_TEST(v[1]->s == "one point zero");
    BOOST_TEST(v[1]->hits == 1);
    BOOST_TEST(dists_sq[1] == 1.0);
}

BOOST_AUTO_TEST_CASE(over_capacity_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    Str s_one("one point zero");
    kb.insert(1.0, &s_one);
    Str s_three("three point zero");
    kb.insert(3.0, &s_three);
    Str s_two("two point zero");
    kb.insert(2.0, &s_two);
    Str s_first_extra("four point four");
    kb.insert(4.4, &s_first_extra);
    Str s_second_extra("two point two");
    kb.insert(2.2, &s_second_extra);
    BOOST_TEST(kb.kth_dist_sq == 2.2);
    BOOST_TEST(kb.hits == 3);
    vector<Str*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_REQUIRE(v.size() == 3);
    BOOST_TEST(dists_sq.size() == 3);
    BOOST_TEST(v[0]->s == "two point two");
    BOOST_TEST(v[0]->hits == 1);
    BOOST_TEST(dists_sq[0] == 2.2);
    BOOST_TEST(v[1]->s == "two point zero");
    BOOST_TEST(v[1]->hits == 1);
    BOOST_TEST(dists_sq[1] == 2.0);
    BOOST_TEST(v[2]->s == "one point zero");
    BOOST_TEST(v[2]->hits == 1);
    BOOST_TEST(dists_sq[2] == 1.0);
}

BOOST_AUTO_TEST_CASE(over_capacity_hits_gt1_test, *tolerance(TOL)) {
    KBest<Str> kb(3);
    Str s_one("one point zero", 2);
    kb.insert(1.0, &s_one);
    Str s_three("three point zero", 2);
    kb.insert(3.0, &s_three);
    Str s_two("two point zero", 2);
    kb.insert(2.0, &s_two);
    BOOST_TEST(kb.kth_dist_sq == 2.0);
    BOOST_TEST(kb.hits == 4);
    vector<Str*> v;
    vector<double> dists_sq;
    kb.fill(&v, &dists_sq);
    BOOST_REQUIRE(v.size() == 2);
    BOOST_TEST(v[0]->s == "two point zero");
    BOOST_TEST(v[0]->hits == 2);
    BOOST_TEST(dists_sq[0] == 2.0);
    BOOST_TEST(v[1]->s == "one point zero");
    BOOST_TEST(v[1]->hits == 2);
    BOOST_TEST(dists_sq[1] == 1.0);
}

BOOST_AUTO_TEST_SUITE_END();  // k_best_suite
BOOST_AUTO_TEST_SUITE_END();  // kd_tree_suite

}  // namespace kd_tree
}  // namespace fluoroseq
