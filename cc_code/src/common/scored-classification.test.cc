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
#include "scored-classification.h"

// Standard C++ library headers:
#include <climits>

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(common_suite)
BOOST_AUTO_TEST_SUITE(scored_classification_suite)

BOOST_AUTO_TEST_CASE(constructor_test) {
    int id = 42;
    double score = 0.9;
    double total = 3.14;
    ScoredClassification sc(id, score, total);
    BOOST_TEST(sc.id == id);
    BOOST_TEST(sc.score == score);
    BOOST_TEST(sc.total == total);
}

BOOST_AUTO_TEST_CASE(constructor_default_test) {
    ScoredClassification sc;
    BOOST_TEST(sc.id == -1);
    BOOST_TEST(sc.score == INT_MIN);
    BOOST_TEST(sc.total == 0.0);
}

BOOST_AUTO_TEST_CASE(adjusted_score_test, *tolerance(TOL)) {
    int id = 42;
    double score = 0.9;
    double total = 3.14;
    ScoredClassification sc(id, score, total);
    BOOST_TEST(sc.adjusted_score() == score / total);
}

BOOST_AUTO_TEST_CASE(gt_op_true_test) {
    int id1 = 421;
    double score1 = 0.9;
    double total1 = 3.14;
    ScoredClassification sc1(id1, score1, total1);
    int id2 = 422;
    double score2 = 0.8;
    double total2 = 3.14;
    ScoredClassification sc2(id2, score2, total2);
    BOOST_TEST((sc1 > sc2));
}

BOOST_AUTO_TEST_CASE(gt_op_false_test) {
    int id1 = 421;
    double score1 = 0.1;
    double total1 = 3.14;
    ScoredClassification sc1(id1, score1, total1);
    int id2 = 422;
    double score2 = 0.2;
    double total2 = 3.14;
    ScoredClassification sc2(id2, score2, total2);
    BOOST_TEST(!(sc1 > sc2));
}

BOOST_AUTO_TEST_SUITE_END()  // scored_classification_suite
BOOST_AUTO_TEST_SUITE_END()  // common_suite

}  // namespace whatprot
