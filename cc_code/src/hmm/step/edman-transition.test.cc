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
#include <cmath>

// File under test:
#include "edman-transition.h"

// Local project headers:
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(edman_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    BOOST_TEST(et.p_edman_failure == p_fail);
    // Should also test that the DyeSeq and DyeTrack were copied over, but this
    // would require equality operators for those classes which I don't want to
    // write right now. Anyways this should be covered by the forward tests.
}

BOOST_AUTO_TEST_CASE(prune_forward_test) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    KDRange range;
    range.min = {1u, 2u, 3u};
    range.max = {3u, 4u, 5u};
    bool allow_detached;
    et.prune_forward(&range, &allow_detached);
    BOOST_TEST(et.true_forward_range.min[0] == 1u);
    BOOST_TEST(et.true_forward_range.min[1] == 2u);
    BOOST_TEST(et.true_forward_range.min[2] == 3u);
    BOOST_TEST(et.true_forward_range.max[0] == 3u);
    BOOST_TEST(et.true_forward_range.max[1] == 4u);
    BOOST_TEST(et.true_forward_range.max[2] == 5u);
    BOOST_TEST(range.min[0] == 1u);
    BOOST_TEST(range.min[1] == 1u);
    BOOST_TEST(range.min[2] == 2u);
    BOOST_TEST(range.max[0] == 4u);
    BOOST_TEST(range.max[1] == 4u);
    BOOST_TEST(range.max[2] == 5u);
    BOOST_TEST(et.safe_backward_range.min[0] == 1u);
    BOOST_TEST(et.safe_backward_range.min[1] == 1u);
    BOOST_TEST(et.safe_backward_range.min[2] == 2u);
    BOOST_TEST(et.safe_backward_range.max[0] == 4u);
    BOOST_TEST(et.safe_backward_range.max[1] == 4u);
    BOOST_TEST(et.safe_backward_range.max[2] == 5u);
}

BOOST_AUTO_TEST_CASE(prune_backward_test) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    KDRange range;
    range.min = {1u, 2u, 3u};
    range.max = {5u, 6u, 7u};
    et.true_forward_range.min = {1u, 3u, 5u};
    et.true_forward_range.max = {4u, 8u, 8u};
    et.safe_backward_range.min = {1u, 2u, 4u};
    et.safe_backward_range.max = {5u, 5u, 7u};
    bool allow_detached;
    et.prune_backward(&range, &allow_detached);
    BOOST_TEST(et.true_backward_range.min[0] == 1u);
    BOOST_TEST(et.true_backward_range.min[1] == 2u);
    BOOST_TEST(et.true_backward_range.min[2] == 4u);
    BOOST_TEST(et.true_backward_range.max[0] == 5u);
    BOOST_TEST(et.true_backward_range.max[1] == 5u);
    BOOST_TEST(et.true_backward_range.max[2] == 7u);
    BOOST_TEST(range.min[0] == 1u);
    BOOST_TEST(range.min[1] == 3u);
    BOOST_TEST(range.min[2] == 5u);
    BOOST_TEST(range.max[0] == 4u);
    BOOST_TEST(range.max[1] == 6u);
    BOOST_TEST(range.max[2] == 8u);
    BOOST_TEST(et.safe_forward_range.min[0] == 0u);
    BOOST_TEST(et.safe_forward_range.min[1] == 2u);
    BOOST_TEST(et.safe_forward_range.min[2] == 4u);
    BOOST_TEST(et.safe_forward_range.max[0] == 5u);
    BOOST_TEST(et.safe_forward_range.max[1] == 6u);
    BOOST_TEST(et.safe_forward_range.max[2] == 8u);
    BOOST_TEST(et.true_forward_range.min[0] == 1u);
    BOOST_TEST(et.true_forward_range.min[1] == 3u);
    BOOST_TEST(et.true_forward_range.min[2] == 5u);
    BOOST_TEST(et.true_forward_range.max[0] == 4u);
    BOOST_TEST(et.true_forward_range.max[1] == 6u);
    BOOST_TEST(et.true_forward_range.max[2] == 8u);
}

BOOST_AUTO_TEST_CASE(forward_trivial_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 1};
    et.safe_forward_range.max = {1, 1};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 1};
    et.safe_backward_range.max = {2, 1};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 1.0;
    psv1.tensor[{1, 0}] = -1000.0;  // to be ignored
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 1.0 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0}]) == 1.0 * p_pop);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_basic_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {1, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.tensor[{1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.7 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0}]) == 0.3 * p_pop);
    BOOST_TEST((psv2->tensor[{1, 1}]) == 0.7 * p_pop);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_more_edmans_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {3, 1};
    et.safe_forward_range.max = {3, 1};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {4, 1};
    et.safe_backward_range.max = {4, 1};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 4;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.2;
    psv1.tensor[{1, 0}] = 0.3;
    psv1.tensor[{2, 0}] = 0.5;
    psv1.tensor[{3, 0}] = -1000.0;  // to be ignored
    unsigned int edmans = 2;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 3u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.2 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0}]) == 0.2 * p_pop + 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{2, 0}]) == 0.3 * p_pop + 0.5 * p_fail);
    BOOST_TEST((psv2->tensor[{3, 0}]) == 0.5 * p_pop);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_dye_colors_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0, 0};
    et.safe_forward_range.min = {0, 0, 0};
    et.true_forward_range.max = {1, 2, 2};
    et.safe_forward_range.max = {1, 2, 2};
    et.true_backward_range.min = {0, 0, 0};
    et.safe_backward_range.min = {0, 0, 0};
    et.true_backward_range.max = {2, 2, 2};
    et.safe_backward_range.max = {2, 2, 2};
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 1}] = 0.2;
    psv1.tensor[{0, 1, 0}] = 0.3;
    psv1.tensor[{0, 1, 1}] = 0.4;
    psv1.tensor[{1, 0, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 0, 1}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0, 0}]) == 0.1 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 0, 1}]) == 0.2 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1, 0}]) == 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1, 1}]) == 0.4 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0, 0}]) == 0.1 * p_pop);
    BOOST_TEST((psv2->tensor[{1, 0, 1}]) == 0.2 * p_pop);
    BOOST_TEST((psv2->tensor[{1, 1, 0}]) == 0.3 * p_pop);
    BOOST_TEST((psv2->tensor[{1, 1, 1}]) == 0.4 * p_pop);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_irrelevant_dye_seq_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, ".0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {1, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.tensor[{1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.7 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0}]) == 0.3 * p_pop);
    BOOST_TEST((psv2->tensor[{1, 1}]) == 0.7 * p_pop);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_one_dye_first_edman_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {1, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.tensor[{1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.7 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0}]) == (0.3 + 0.7) * p_pop);
    BOOST_TEST((psv2->tensor[{1, 1}]) == 0.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_two_dyes_second_edman_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "00");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {2, 3};
    et.safe_forward_range.max = {2, 3};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {3, 3};
    et.safe_backward_range.max = {3, 3};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.1;
    psv1.tensor[{0, 1}] = 0.2;
    psv1.tensor[{0, 2}] = 0.3;
    psv1.tensor[{1, 0}] = 0.4;
    psv1.tensor[{1, 1}] = 0.5;
    psv1.tensor[{1, 2}] = 0.0;
    psv1.tensor[{2, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{2, 1}] = -1000.0;  // to be ignored
    psv1.tensor[{2, 2}] = -1000.0;  // to be ignored
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 2u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.1 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.2 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 2}]) == 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0}])
               == (0.1 + 0.2 / 2.0) * p_pop + 0.4 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 1}])
               == (0.2 / 2.0 + 0.3) * p_pop + 0.5 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 2}]) == 0.0);
    BOOST_TEST((psv2->tensor[{2, 0}]) == (0.4 + 0.5) * p_pop);
    BOOST_TEST((psv2->tensor[{2, 1}]) == 0.0);
    BOOST_TEST((psv2->tensor[{2, 2}]) == 0.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_three_dyes_first_edman_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "000");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 4};
    et.safe_forward_range.max = {1, 4};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 4};
    et.safe_backward_range.max = {2, 4};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 4;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.1;
    psv1.tensor[{0, 1}] = 0.2;
    psv1.tensor[{0, 2}] = 0.3;
    psv1.tensor[{0, 3}] = 0.4;
    psv1.tensor[{1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 2}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 3}] = -1000.0;  // to be ignored
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.1 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.2 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 2}]) == 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 3}]) == 0.4 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0}]) == (0.1 + 0.2 / 3.0) * p_pop);
    BOOST_TEST((psv2->tensor[{1, 1}])
               == (0.2 * 2.0 / 3.0 + 0.3 * 2.0 / 3.0) * p_pop);
    BOOST_TEST((psv2->tensor[{1, 2}]) == (0.3 / 3.0 + 0.4) * p_pop);
    BOOST_TEST((psv2->tensor[{1, 3}]) == 0.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_two_dye_colors_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0, 0};
    et.safe_forward_range.min = {0, 0, 0};
    et.true_forward_range.max = {2, 2, 2};
    et.safe_forward_range.max = {2, 2, 2};
    et.true_backward_range.min = {0, 0, 0};
    et.safe_backward_range.min = {0, 0, 0};
    et.true_backward_range.max = {3, 2, 2};
    et.safe_backward_range.max = {3, 2, 2};
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 1}] = 0.2;
    psv1.tensor[{0, 1, 0}] = 0.3;
    psv1.tensor[{0, 1, 1}] = 0.4;
    psv1.tensor[{1, 0, 0}] = 0.5;
    psv1.tensor[{1, 0, 1}] = 0.6;
    psv1.tensor[{1, 1, 0}] = 0.0;
    psv1.tensor[{1, 1, 1}] = 0.0;
    psv1.tensor[{2, 0, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{2, 0, 1}] = -1000.0;  // to be ignored
    psv1.tensor[{2, 1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{2, 1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 2u);
    BOOST_TEST((psv2->tensor[{0, 0, 0}]) == 0.1 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 0, 1}]) == 0.2 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1, 0}]) == 0.3 * p_fail);
    BOOST_TEST((psv2->tensor[{0, 1, 1}]) == 0.4 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0, 0}]) == (0.1 + 0.3) * p_pop + 0.5 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 0, 1}]) == (0.2 + 0.4) * p_pop + 0.6 * p_fail);
    BOOST_TEST((psv2->tensor[{1, 1, 0}]) == 0.0);
    BOOST_TEST((psv2->tensor[{1, 1, 1}]) == 0.0);
    BOOST_TEST((psv2->tensor[{2, 0, 0}]) == (0.5 + 0.6) * p_pop);
    BOOST_TEST((psv2->tensor[{2, 0, 1}]) == 0.0);
    BOOST_TEST((psv2->tensor[{2, 1, 0}]) == 0.0);
    BOOST_TEST((psv2->tensor[{2, 1, 1}]) == 0.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_with_empty_forward_range_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {0, 0};
    et.safe_forward_range.max = {0, 0};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{0, 1}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.0);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.0);
    BOOST_TEST((psv2->tensor[{1, 0}]) == 0.0);
    BOOST_TEST((psv2->tensor[{1, 1}]) == 0.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_no_crash_with_uninitialized_out_of_range,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {1, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {0, 0};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 1.00;
    psv1.tensor[{0, 1}] = 1.01;
    psv1.tensor[{1, 0}] = 1.10;
    psv1.tensor[{1, 1}] = 1.11;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = et.forward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_trivial_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 1};
    et.safe_forward_range.max = {2, 1};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 1};
    et.safe_backward_range.max = {2, 1};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{1, 0}] = 0.7;
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_fail * 0.3 + p_pop * 0.7);
    // We ignore the value at {1, 0}, it doesn't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_basic_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {2, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.tensor[{1, 0}] = 1.33;
    psv1.tensor[{1, 1}] = 1.77;
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_fail * 0.3 + p_pop * 1.33);
    BOOST_TEST((psv2->tensor[{0, 1}]) == p_fail * 0.7 + p_pop * 1.77);
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_more_edmans_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {3, 1};
    et.safe_forward_range.max = {4, 1};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {4, 1};
    et.safe_backward_range.max = {4, 1};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 4;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.2;
    psv1.tensor[{1, 0}] = 0.3;
    psv1.tensor[{2, 0}] = 0.5;
    psv1.tensor[{3, 0}] = 0.7;
    unsigned int edmans = 3;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 2u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_fail * 0.2 + p_pop * 0.3);
    BOOST_TEST((psv2->tensor[{1, 0}]) == p_fail * 0.3 + p_pop * 0.5);
    BOOST_TEST((psv2->tensor[{2, 0}]) == p_fail * 0.5 + p_pop * 0.7);
    // We ignore the value at {3, 0}, it doesn't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_multiple_dye_colors_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0, 0};
    et.safe_forward_range.min = {0, 0, 0};
    et.true_forward_range.max = {1, 2, 2};
    et.safe_forward_range.max = {2, 2, 2};
    et.true_backward_range.min = {0, 0, 0};
    et.safe_backward_range.min = {0, 0, 0};
    et.true_backward_range.max = {2, 2, 2};
    et.safe_backward_range.max = {2, 2, 2};
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 1}] = 0.2;
    psv1.tensor[{0, 1, 0}] = 0.3;
    psv1.tensor[{0, 1, 1}] = 0.4;
    psv1.tensor[{1, 0, 0}] = 1.11;
    psv1.tensor[{1, 0, 1}] = 1.22;
    psv1.tensor[{1, 1, 0}] = 1.33;
    psv1.tensor[{1, 1, 1}] = 1.44;
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    BOOST_TEST((psv2->tensor[{0, 0, 0}]) == p_fail * 0.1 + p_pop * 1.11);
    BOOST_TEST((psv2->tensor[{0, 0, 1}]) == p_fail * 0.2 + p_pop * 1.22);
    BOOST_TEST((psv2->tensor[{0, 1, 0}]) == p_fail * 0.3 + p_pop * 1.33);
    BOOST_TEST((psv2->tensor[{0, 1, 1}]) == p_fail * 0.4 + p_pop * 1.44);
    // We ignore the values at {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, and {1, 1, 1};
    // they don't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_irrelevant_dye_seq_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, ".0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {2, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.tensor[{1, 0}] = 1.33;
    psv1.tensor[{1, 1}] = 1.77;
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_fail * 0.3 + p_pop * 1.33);
    BOOST_TEST((psv2->tensor[{0, 1}]) == p_fail * 0.7 + p_pop * 1.77);
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_one_dye_first_edman_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {2, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.tensor[{1, 0}] = 1.33;
    psv1.tensor[{1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_fail * 0.3 + p_pop * 1.33);
    BOOST_TEST((psv2->tensor[{0, 1}]) == p_fail * 0.7 + p_pop * 1.33);
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_two_dyes_second_edman_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "00");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {2, 3};
    et.safe_forward_range.max = {3, 3};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {3, 3};
    et.safe_backward_range.max = {3, 3};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.1;
    psv1.tensor[{0, 1}] = 0.2;
    psv1.tensor[{0, 2}] = 0.3;
    psv1.tensor[{1, 0}] = 0.4;
    psv1.tensor[{1, 1}] = 0.5;
    psv1.tensor[{1, 2}] = 0.0;
    psv1.tensor[{2, 0}] = 0.6;
    psv1.tensor[{2, 1}] = 0.0;
    psv1.tensor[{2, 2}] = 0.0;
    unsigned int edmans = 2;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_fail * 0.1 + p_pop * 0.4);
    BOOST_TEST((psv2->tensor[{0, 1}])
               == p_fail * 0.2 + p_pop * (0.5 / 2.0 + 0.4 / 2.0));
    BOOST_TEST((psv2->tensor[{0, 2}]) == p_fail * 0.3 + p_pop * 0.5);
    BOOST_TEST((psv2->tensor[{1, 0}]) == p_fail * 0.4 + p_pop * 0.6);
    BOOST_TEST((psv2->tensor[{1, 1}]) == p_fail * 0.5 + p_pop * 0.6);
    BOOST_TEST((psv2->tensor[{1, 2}]) == 0.0);
    // We ignore the values at {2, 0}, {2, 1}, and {2, 2}. They don't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_three_dyes_first_edman_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "000");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 4};
    et.safe_forward_range.max = {2, 4};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 4};
    et.safe_backward_range.max = {2, 4};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 4;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.1;
    psv1.tensor[{0, 1}] = 0.2;
    psv1.tensor[{0, 2}] = 0.3;
    psv1.tensor[{0, 3}] = 0.4;
    psv1.tensor[{1, 0}] = 1.11;
    psv1.tensor[{1, 1}] = 1.22;
    psv1.tensor[{1, 2}] = 1.33;
    psv1.tensor[{1, 3}] = -1000.0;  // to be ignored
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == p_fail * 0.1 + p_pop * 1.11);
    BOOST_TEST((psv2->tensor[{0, 1}])
               == p_fail * 0.2 + p_pop * (1.11 * 1.0 / 3.0 + 1.22 * 2.0 / 3.0));
    BOOST_TEST((psv2->tensor[{0, 2}])
               == p_fail * 0.3 + p_pop * (1.22 * 2.0 / 3.0 + 1.33 * 1.0 / 3.0));
    BOOST_TEST((psv2->tensor[{0, 3}]) == p_fail * 0.4 + p_pop * 1.33);
    // We ignore the values at {1, 0}, {1, 1}, {1, 2}, and {1, 3}.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_two_dye_colors_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0, 0};
    et.safe_forward_range.min = {0, 0, 0};
    et.true_forward_range.max = {2, 2, 2};
    et.safe_forward_range.max = {3, 2, 2};
    et.true_backward_range.min = {0, 0, 0};
    et.safe_backward_range.min = {0, 0, 0};
    et.true_backward_range.max = {3, 2, 2};
    et.safe_backward_range.max = {3, 2, 2};
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 1}] = 0.2;
    psv1.tensor[{0, 1, 0}] = 0.3;
    psv1.tensor[{0, 1, 1}] = 0.4;
    psv1.tensor[{1, 0, 0}] = 0.5;
    psv1.tensor[{1, 0, 1}] = 0.6;
    psv1.tensor[{1, 1, 0}] = 0.0;
    psv1.tensor[{1, 1, 1}] = 0.0;
    psv1.tensor[{2, 0, 0}] = 0.7;
    psv1.tensor[{2, 0, 1}] = 0.0;
    psv1.tensor[{2, 1, 0}] = 0.0;
    psv1.tensor[{2, 1, 1}] = 0.0;
    unsigned int edmans = 2;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 1u);
    BOOST_TEST((psv2->tensor[{0, 0, 0}]) == p_fail * 0.1 + p_pop * 0.5);
    BOOST_TEST((psv2->tensor[{0, 0, 1}]) == p_fail * 0.2 + p_pop * 0.6);
    BOOST_TEST((psv2->tensor[{0, 1, 0}]) == p_fail * 0.3 + p_pop * 0.5);
    BOOST_TEST((psv2->tensor[{0, 1, 1}]) == p_fail * 0.4 + p_pop * 0.6);
    BOOST_TEST((psv2->tensor[{1, 0, 0}]) == p_fail * 0.5 + p_pop * 0.7);
    BOOST_TEST((psv2->tensor[{1, 0, 1}]) == p_fail * 0.6 + p_pop * 0.7);
    BOOST_TEST((psv2->tensor[{1, 1, 0}]) == 0.0);
    BOOST_TEST((psv2->tensor[{1, 1, 1}]) == 0.0);
    // We ignore the values at {2, 0, 0}, {2, 0, 1}, {2, 1, 0}, and {2, 1, 1}.
    // They don't matter.
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_with_empty_backward_range_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {1, 2};
    et.safe_forward_range.max = {1, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {0, 0};
    et.safe_backward_range.max = {0, 0};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{0, 1}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 0}] = -1000.0;  // to be ignored
    psv1.tensor[{1, 1}] = -1000.0;  // to be ignored
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 0.0);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 0.0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_no_crash_with_uninitialized_out_of_range,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    et.true_forward_range.min = {0, 0};
    et.safe_forward_range.min = {0, 0};
    et.true_forward_range.max = {0, 0};
    et.safe_forward_range.max = {2, 2};
    et.true_backward_range.min = {0, 0};
    et.safe_backward_range.min = {0, 0};
    et.true_backward_range.max = {2, 2};
    et.safe_backward_range.max = {2, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 1.00;
    psv1.tensor[{0, 1}] = 1.01;
    psv1.tensor[{1, 0}] = 1.10;
    psv1.tensor[{1, 1}] = 1.11;
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = et.backward(psv1, &edmans);
    BOOST_TEST(edmans == 0u);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.61;
    fpsv.tensor[{0, 1}] = 0.91;
    fpsv.tensor[{1, 0}] = 0.62;
    fpsv.tensor[{1, 1}] = 0.92;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.51;
    bpsv.tensor[{0, 1}] = 0.81;
    bpsv.tensor[{1, 0}] = 0.52;
    bpsv.tensor[{1, 1}] = 0.82;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.41;
    nbpsv.tensor[{0, 1}] = 0.71;
    nbpsv.tensor[{1, 0}] = 0.42;
    nbpsv.tensor[{1, 1}] = 0.72;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 1.0;
    SequencingModelFitter smf;
    et.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    BOOST_TEST(smf.p_edman_failure_fit.get()
               == (0.91 * p_fail * 0.71) / (0.91 * 0.81));
}

BOOST_AUTO_TEST_CASE(improve_fit_twice_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector fpsv1(order, shape);
    fpsv1.tensor[{0, 0}] = 0.31;
    fpsv1.tensor[{0, 1}] = 0.91;
    fpsv1.tensor[{1, 0}] = 0.32;
    fpsv1.tensor[{1, 1}] = 0.92;
    PeptideStateVector bpsv1(order, shape);
    bpsv1.tensor[{0, 0}] = 0.331;
    bpsv1.tensor[{0, 1}] = 0.81;
    bpsv1.tensor[{1, 0}] = 0.332;
    bpsv1.tensor[{1, 1}] = 0.82;
    PeptideStateVector nbpsv1(order, shape);
    nbpsv1.tensor[{0, 0}] = 0.21;
    nbpsv1.tensor[{0, 1}] = 0.71;
    nbpsv1.tensor[{1, 0}] = 0.22;
    nbpsv1.tensor[{1, 1}] = 0.72;
    PeptideStateVector fpsv2(order, shape);
    fpsv2.tensor[{0, 0}] = 0.221;
    fpsv2.tensor[{0, 1}] = 0.61;
    fpsv2.tensor[{1, 0}] = 0.222;
    fpsv2.tensor[{1, 1}] = 0.92;
    PeptideStateVector bpsv2(order, shape);
    bpsv2.tensor[{0, 0}] = 0.11;
    bpsv2.tensor[{0, 1}] = 0.51;
    bpsv2.tensor[{1, 0}] = 0.12;
    bpsv2.tensor[{1, 1}] = 0.82;
    PeptideStateVector nbpsv2(order, shape);
    nbpsv2.tensor[{0, 0}] = 0.111;
    nbpsv2.tensor[{0, 1}] = 0.41;
    nbpsv2.tensor[{1, 0}] = 0.112;
    nbpsv2.tensor[{1, 1}] = 0.72;
    delete[] shape;
    unsigned int edmans = 0;
    double prob1 = 0.12345;
    double prob2 = 0.98765;
    SequencingModelFitter smf;
    et.improve_fit(fpsv1, bpsv1, nbpsv1, edmans, prob1, &smf);
    et.improve_fit(fpsv2, bpsv2, nbpsv2, edmans, prob2, &smf);
    BOOST_TEST(smf.p_edman_failure_fit.get()
               == (0.91 * p_fail * 0.71 / prob1 + 0.61 * p_fail * 0.41 / prob2)
                          / (0.91 * 0.81 / prob1 + 0.61 * 0.51 / prob2));
}

BOOST_AUTO_TEST_SUITE_END()  // edman_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
