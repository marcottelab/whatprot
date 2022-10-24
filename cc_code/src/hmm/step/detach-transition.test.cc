/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// Standard C++ library headers:
#include <vector>

// File under test:
#include "detach-transition.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "util/kd-range.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using std::vector;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(detach_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    BOOST_CHECK_EQUAL(dt.p_detach, p_detach);
}

BOOST_AUTO_TEST_CASE(prune_forward_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    KDRange range;
    range.min = vector<unsigned int>(2, 1u);
    range.max = vector<unsigned int>(2, 3u);
    bool allow_detached;
    dt.prune_forward(&range, &allow_detached);
    BOOST_TEST(dt.pruned_range.min[0] == 1u);
    BOOST_TEST(dt.pruned_range.min[1] == 1u);
    BOOST_TEST(dt.pruned_range.max[0] == 3u);
    BOOST_TEST(dt.pruned_range.max[1] == 3u);
    BOOST_TEST(range.min[0] == 1u);
    BOOST_TEST(range.min[1] == 1u);
    BOOST_TEST(range.max[0] == 3u);
    BOOST_TEST(range.max[1] == 3u);
}

BOOST_AUTO_TEST_CASE(prune_backward_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    KDRange range;
    range.min = vector<unsigned int>(2, 1u);
    range.max = vector<unsigned int>(2, 3u);
    dt.pruned_range.min = vector<unsigned int>(2, 2u);
    dt.pruned_range.max = vector<unsigned int>(2, 4u);
    bool allow_detached;
    dt.prune_backward(&range, &allow_detached);
    BOOST_TEST(dt.pruned_range.min[0] == 2u);
    BOOST_TEST(dt.pruned_range.min[1] == 2u);
    BOOST_TEST(dt.pruned_range.max[0] == 3u);
    BOOST_TEST(dt.pruned_range.max[1] == 3u);
    BOOST_TEST(range.min[0] == 2u);
    BOOST_TEST(range.min[1] == 2u);
    BOOST_TEST(range.max[0] == 3u);
    BOOST_TEST(range.max[1] == 3u);
}

BOOST_AUTO_TEST_CASE(forward_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 1};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 1.0;
    psv1.broken_n_tensor[{0, 0}] = 2.0;
    psv1.p_detached = 1.0;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 1.0 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 2.0 * (1 - p_detach));
    BOOST_TEST(psv2->p_detached == 1.0 + (1.0 + 2.0) * p_detach);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.07;
    psv1.p_detached = 0.9;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.7 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.03 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.07 * (1 - p_detach));
    BOOST_TEST(psv2->p_detached == (0.3 + 0.7 + 0.03 + 0.07) * p_detach + 0.9);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.6 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 2}]) == 0.1 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.03 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.06 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}]) == 0.01 * (1 - p_detach));
    BOOST_TEST(psv2->p_detached
               == (0.3 + 0.6 + 0.1 + 0.03 + 0.06 + 0.01) * p_detach + 0.2);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {3, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.1;
    psv1.tensor[{0, 1}] = 0.2;
    psv1.tensor[{1, 0}] = 0.3;
    psv1.tensor[{1, 1}] = 0.4;
    psv1.tensor[{2, 0}] = 0.5;
    psv1.tensor[{2, 1}] = 0.6;
    psv1.broken_n_tensor[{0, 0}] = 0.01;
    psv1.broken_n_tensor[{0, 1}] = 0.02;
    psv1.broken_n_tensor[{1, 0}] = 0.03;
    psv1.broken_n_tensor[{1, 1}] = 0.04;
    psv1.broken_n_tensor[{2, 0}] = 0.05;
    psv1.broken_n_tensor[{2, 1}] = 0.06;
    psv1.p_detached = 0.7;
    unsigned int edmans = 2;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.2 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{1, 1}]) == 0.4 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{2, 1}]) == 0.6 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.02 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{1, 1}]) == 0.04 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{2, 1}]) == 0.06 * (1 - p_detach));
    // Distribution between empty states is of no importance. We only care about
    // the sum.
    double sum_empties = 0.0;
    sum_empties += psv2->tensor[{0, 0}];
    sum_empties += psv2->tensor[{1, 0}];
    sum_empties += psv2->tensor[{2, 0}];
    sum_empties += psv2->broken_n_tensor[{0, 0}];
    sum_empties += psv2->broken_n_tensor[{1, 0}];
    sum_empties += psv2->broken_n_tensor[{2, 0}];
    sum_empties += psv2->p_detached;
    BOOST_TEST(sum_empties
               == 0.1 + 0.3 + 0.5 + 0.01 + 0.03 + 0.05
                          + (0.2 + 0.4 + 0.6 + 0.02 + 0.04 + 0.06) * p_detach
                          + 0.7);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_dye_colors_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0, 0};
    dt.pruned_range.max = {1, 2, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 1}] = 0.2;
    psv1.tensor[{0, 1, 0}] = 0.3;
    psv1.tensor[{0, 1, 1}] = 0.4;
    psv1.broken_n_tensor[{0, 0, 0}] = 0.01;
    psv1.broken_n_tensor[{0, 0, 1}] = 0.02;
    psv1.broken_n_tensor[{0, 1, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1, 1}] = 0.04;
    psv1.p_detached = 0.5;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0, 0}]) == 0.1 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 0, 1}]) == 0.2 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 1, 0}]) == 0.3 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 1, 1}]) == 0.4 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0}]) == 0.01 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 1}]) == 0.02 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 0}]) == 0.03 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 1}]) == 0.04 * (1 - p_detach));
    BOOST_TEST(psv2->p_detached
               == (0.1 + 0.2 + 0.3 + 0.4 + 0.01 + 0.02 + 0.03 + 0.04) * p_detach
                          + 0.5);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_pruned_range_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 1};
    dt.pruned_range.max = {1, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.6 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.06 * (1 - p_detach));
    BOOST_TEST(psv2->p_detached == (0.6 + 0.06) * p_detach + 0.2);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_no_detached_forward_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = false;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.6 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 2}]) == 0.1 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.03 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.06 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}]) == 0.01 * (1 - p_detach));
    BOOST_TEST(psv2->p_detached
               == (0.3 + 0.6 + 0.1 + 0.03 + 0.06 + 0.01) * p_detach);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_no_detached_backward_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = true;
    dt.detached_backward = false;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.6 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 2}]) == 0.1 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.03 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.06 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}]) == 0.01 * (1 - p_detach));
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_no_detached_forward_or_backward_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = false;
    dt.detached_backward = false;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.6 * (1 - p_detach));
    BOOST_TEST((psv2->tensor[{0, 2}]) == 0.1 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.03 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.06 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}]) == 0.01 * (1 - p_detach));
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 1};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 1.0;
    psv1.broken_n_tensor[{0, 0}] = 2.0;
    psv1.p_detached = 1.0;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 1.0 * p_detach + 1.0 * (1 - p_detach));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == 1.0 * p_detach + 2.0 * (1 - p_detach));
    BOOST_TEST(psv2->p_detached == 1.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.07;
    psv1.p_detached = 0.9;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_detach * 0.9 + (1 - p_detach) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1}]) == p_detach * 0.9 + (1 - p_detach) * 0.7);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == p_detach * 0.9 + (1 - p_detach) * 0.03);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == p_detach * 0.9 + (1 - p_detach) * 0.07);
    BOOST_TEST(psv2->p_detached = 0.9);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_detach * 0.2 + (1 - p_detach) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1}]) == p_detach * 0.2 + (1 - p_detach) * 0.6);
    BOOST_TEST((psv2->tensor[{0, 2}]) == p_detach * 0.2 + (1 - p_detach) * 0.1);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == p_detach * 0.2 + (1 - p_detach) * 0.03);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == p_detach * 0.2 + (1 - p_detach) * 0.06);
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}])
               == p_detach * 0.2 + (1 - p_detach) * 0.01);
    BOOST_TEST(psv2->p_detached = 0.2);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {3, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.88;
    psv1.tensor[{0, 1}] = 0.2;
    psv1.tensor[{1, 0}] = 0.88;
    psv1.tensor[{1, 1}] = 0.4;
    psv1.tensor[{2, 0}] = 0.88;
    psv1.tensor[{2, 1}] = 0.6;
    psv1.broken_n_tensor[{0, 0}] = 0.88;
    psv1.broken_n_tensor[{0, 1}] = 0.02;
    psv1.broken_n_tensor[{1, 0}] = 0.88;
    psv1.broken_n_tensor[{1, 1}] = 0.04;
    psv1.broken_n_tensor[{2, 0}] = 0.88;
    psv1.broken_n_tensor[{2, 1}] = 0.06;
    psv1.p_detached = 0.88;
    unsigned int edmans = 2;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.88);
    BOOST_TEST((psv2->tensor[{0, 1}])
               == p_detach * 0.88 + (1 - p_detach) * 0.2);
    BOOST_TEST((psv2->tensor[{1, 0}]) == 0.88);
    BOOST_TEST((psv2->tensor[{1, 1}])
               == p_detach * 0.88 + (1 - p_detach) * 0.4);
    BOOST_TEST((psv2->tensor[{2, 0}]) == 0.88);
    BOOST_TEST((psv2->tensor[{2, 1}])
               == p_detach * 0.88 + (1 - p_detach) * 0.6);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.88);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == p_detach * 0.88 + (1 - p_detach) * 0.02);
    BOOST_TEST((psv2->broken_n_tensor[{1, 0}]) == 0.88);
    BOOST_TEST((psv2->broken_n_tensor[{1, 1}])
               == p_detach * 0.88 + (1 - p_detach) * 0.04);
    BOOST_TEST((psv2->broken_n_tensor[{2, 0}]) == 0.88);
    BOOST_TEST((psv2->broken_n_tensor[{2, 1}])
               == p_detach * 0.88 + (1 - p_detach) * 0.06);
    BOOST_TEST(psv2->p_detached == 0.88);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_multiple_dye_colors_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0, 0};
    dt.pruned_range.max = {1, 2, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 1}] = 0.2;
    psv1.tensor[{0, 1, 0}] = 0.3;
    psv1.tensor[{0, 1, 1}] = 0.4;
    psv1.broken_n_tensor[{0, 0, 0}] = 0.01;
    psv1.broken_n_tensor[{0, 0, 1}] = 0.02;
    psv1.broken_n_tensor[{0, 1, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1, 1}] = 0.04;
    psv1.p_detached = 0.5;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0, 0}])
               == p_detach * 0.5 + (1 - p_detach) * 0.1);
    BOOST_TEST((psv2->tensor[{0, 0, 1}])
               == p_detach * 0.5 + (1 - p_detach) * 0.2);
    BOOST_TEST((psv2->tensor[{0, 1, 0}])
               == p_detach * 0.5 + (1 - p_detach) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1, 1}])
               == p_detach * 0.5 + (1 - p_detach) * 0.4);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0}])
               == p_detach * 0.5 + (1 - p_detach) * 0.01);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 1}])
               == p_detach * 0.5 + (1 - p_detach) * 0.02);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 0}])
               == p_detach * 0.5 + (1 - p_detach) * 0.03);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 1}])
               == p_detach * 0.5 + (1 - p_detach) * 0.04);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_pruned_range_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 1};
    dt.pruned_range.max = {1, 2};
    dt.detached_forward = true;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 1}]) == p_detach * 0.2 + (1 - p_detach) * 0.6);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == p_detach * 0.2 + (1 - p_detach) * 0.06);
    BOOST_TEST(psv2->p_detached == 0.2);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_no_detached_forward_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = false;
    dt.detached_backward = true;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_detach * 0.2 + (1 - p_detach) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1}]) == p_detach * 0.2 + (1 - p_detach) * 0.6);
    BOOST_TEST((psv2->tensor[{0, 2}]) == p_detach * 0.2 + (1 - p_detach) * 0.1);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == p_detach * 0.2 + (1 - p_detach) * 0.03);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == p_detach * 0.2 + (1 - p_detach) * 0.06);
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}])
               == p_detach * 0.2 + (1 - p_detach) * 0.01);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_no_detached_backward_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = true;
    dt.detached_backward = false;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == (1 - p_detach) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1}]) == (1 - p_detach) * 0.6);
    BOOST_TEST((psv2->tensor[{0, 2}]) == (1 - p_detach) * 0.1);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == (1 - p_detach) * 0.03);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == (1 - p_detach) * 0.06);
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}]) == (1 - p_detach) * 0.01);
    BOOST_TEST(psv2->p_detached == 0.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_no_detached_forward_or_backward_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    dt.detached_forward = false;
    dt.detached_backward = false;
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.6;
    psv1.tensor[{0, 2}] = 0.1;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.06;
    psv1.broken_n_tensor[{0, 2}] = 0.01;
    psv1.p_detached = 0.2;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = dt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == (1 - p_detach) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1}]) == (1 - p_detach) * 0.6);
    BOOST_TEST((psv2->tensor[{0, 2}]) == (1 - p_detach) * 0.1);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == (1 - p_detach) * 0.03);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == (1 - p_detach) * 0.06);
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}]) == (1 - p_detach) * 0.01);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    dt.pruned_range.min = {0, 0};
    dt.pruned_range.max = {1, 3};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    fpsv.tensor[{0, 2}] = 0.91;
    fpsv.broken_n_tensor[{0, 0}] = 0.031;
    fpsv.broken_n_tensor[{0, 1}] = 0.071;
    fpsv.broken_n_tensor[{0, 2}] = 0.091;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    bpsv.tensor[{0, 2}] = 0.92;
    bpsv.broken_n_tensor[{0, 0}] = 0.032;
    bpsv.broken_n_tensor[{0, 1}] = 0.072;
    bpsv.broken_n_tensor[{0, 2}] = 0.092;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    nbpsv.tensor[{0, 2}] = 0.93;
    nbpsv.broken_n_tensor[{0, 0}] = 0.033;
    nbpsv.broken_n_tensor[{0, 1}] = 0.073;
    nbpsv.broken_n_tensor[{0, 2}] = 0.093;
    nbpsv.p_detached = 0.33;  // same as for {0, 0} which will always be true.
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 0.31 * 0.32 + 0.71 * 0.72 + 0.91 * 0.92 + 0.031 * 0.032
                         + 0.071 * 0.072 + 0.091 * 0.092;
    SequencingModelFitter smf;
    dt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    BOOST_TEST(smf.p_detach_fit.get()
               == (0.71 * p_detach * 0.33 + 0.91 * p_detach * 0.33
                   + 0.071 * p_detach * 0.33 + 0.091 * p_detach * 0.33)
                          / (0.71 * 0.72 + 0.91 * 0.92 + 0.071 * 0.072
                             + 0.091 * 0.092));
}

BOOST_AUTO_TEST_SUITE_END()  // detach_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
