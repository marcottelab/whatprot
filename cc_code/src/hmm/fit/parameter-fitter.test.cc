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
#include "parameter-fitter.h"

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(fit_suite)
BOOST_AUTO_TEST_SUITE(parameter_fitter_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    ParameterFitter pf;
    BOOST_TEST(pf.numerator == 0.0);
    BOOST_TEST(pf.denominator == 0.0);
}

BOOST_AUTO_TEST_CASE(get_test, *tolerance(TOL)) {
    ParameterFitter pf;
    pf.numerator = 2.71828;
    pf.denominator = 3.14159;
    BOOST_TEST(pf.get() == 2.71828 / 3.14159);
}

BOOST_AUTO_TEST_SUITE_END()  // parameter_fitter_suite
BOOST_AUTO_TEST_SUITE_END()  // fit_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace fluoroseq
