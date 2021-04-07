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
#include "generic-hmm.h"

// Standard C++ library headers:
#include <cmath>
#include <functional>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/fit/error-model-fitter.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/step/binomial-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/step/edman-transition.h"
#include "tensor/tensor.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using std::exp;
using std::function;
using std::log;
using std::sqrt;
using std::vector;
const double PI = 3.141592653589793238;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(hmm_suite)

BOOST_AUTO_TEST_SUITE_END()  // hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
