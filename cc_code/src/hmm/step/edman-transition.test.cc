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
#include "edman-transition.h"

// Local project headers:
#include "common/error-model.h"
#include "hmm/fit/error-model-fitter.h"

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
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    BOOST_TEST(et.p_edman_failure == p_fail);
    // Should also test that the DyeSeq and DyeTrack were copied over, but this
    // would require equality operators for those classes which I don't want to
    // write right now. Anyways this should be covered by the forward tests.
}

BOOST_AUTO_TEST_CASE(forward_in_place_trivial_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 1.0;  // loc is {0, 0}
    loc[0] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
    int edmans = 0;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 1.0 * p_fail);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(tsr[loc] == 1.0 * p_pop);  // loc is {1, 0}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_trivial_test, *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 1;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, "");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 2;
//     shape[1] = 1;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     tsr1[loc] = 1.0;  // loc is {0, 0}
//     loc[0] = 1;
//     tsr1[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 1.0 * p_fail);  // loc is {0, 0}
//     loc[0] = 1;
//     BOOST_TEST(tsr2[loc] == 1.0 * p_pop);  // loc is {1, 0}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_basic_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
    int edmans = 0;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * p_fail);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * p_pop);  // loc is {1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_basic_test, *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 1;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, "");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 2;
//     shape[1] = 2;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     tsr1[loc] = 0.3;  // loc is {0, 0}
//     loc[1] = 1;
//     tsr1[loc] = 0.7;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     tsr1[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
//     loc[1] = 1;
//     tsr1[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.7 * p_fail);  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.7 * p_pop);  // loc is {1, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_more_edmans_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 3;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 4;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.2;  // loc is {0, 0}
    loc[0] = 1;
    tsr[loc] = 0.3;  // loc is {1, 0}
    loc[0] = 2;
    tsr[loc] = 0.5;  // loc is {2, 0}
    loc[0] = 3;
    tsr[loc] = -1000.0;  // loc is {3, 0} -- this value should be ignored.
    int edmans = 2;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 3);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.2 * p_fail);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * p_pop + 0.3 * p_fail);  // loc is {1, 0}
    loc[0] = 2;
    BOOST_TEST(tsr[loc] == 0.3 * p_pop + 0.5 * p_fail);  // loc is {2, 0}
    loc[0] = 3;
    BOOST_TEST(tsr[loc] == 0.5 * p_pop);  // loc is {3, 0}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_more_edmans_test, *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 3;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, "");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 4;
//     shape[1] = 1;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     tsr1[loc] = 0.2;  // loc is {0, 0}
//     loc[0] = 1;
//     tsr1[loc] = 0.3;  // loc is {1, 0}
//     loc[0] = 2;
//     tsr1[loc] = 0.5;  // loc is {2, 0}
//     loc[0] = 3;
//     tsr1[loc] = -1000.0;  // loc is {3, 0} -- this value should be ignored.
//     int edmans = 2;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 3);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.2 * p_fail);  // loc is {0, 0}
//     loc[0] = 1;
//     BOOST_TEST(tsr2[loc] == 0.2 * p_pop + 0.3 * p_fail);  // loc is {1, 0}
//     loc[0] = 2;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_pop + 0.5 * p_fail);  // loc is {2, 0}
//     loc[0] = 3;
//     BOOST_TEST(tsr2[loc] == 0.5 * p_pop);  // loc is {3, 0}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {1, 0, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 0, 1} -- this value should be ignored.
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {1, 1, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 1, 1} -- this value should be ignored.
    int edmans = 0;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.1 * p_pop);  // loc is {1, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * p_pop);  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_pop);  // loc is {1, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.4 * p_pop);  // loc is {1, 1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_multiple_dye_colors_test,
//                      *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 1;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, "");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 3;
//     int* shape = new int[order];
//     shape[0] = 2;
//     shape[1] = 2;
//     shape[2] = 2;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     tsr1[loc] = 0.1;  // loc is {0, 0, 0}
//     loc[2] = 1;
//     tsr1[loc] = 0.2;  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     tsr1[loc] = 0.3;  // loc is {0, 1, 0}
//     loc[2] = 1;
//     tsr1[loc] = 0.4;  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;
//     tsr1[loc] = -1000.0;  // loc is {1, 0, 0} -- this value should be ignored.
//     loc[2] = 1;
//     tsr1[loc] = -1000.0;  // loc is {1, 0, 1} -- this value should be ignored.
//     loc[1] = 1;
//     loc[2] = 0;
//     tsr1[loc] = -1000.0;  // loc is {1, 1, 0} -- this value should be ignored.
//     loc[2] = 1;
//     tsr1[loc] = -1000.0;  // loc is {1, 1, 1} -- this value should be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.1 * p_pop);  // loc is {1, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.2 * p_pop);  // loc is {1, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_pop);  // loc is {1, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.4 * p_pop);  // loc is {1, 1, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_irrelevant_dye_seq_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, ".0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
    int edmans = 0;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * p_fail);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * p_pop);  // loc is {1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_irrelevant_dye_seq_test, *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 1;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, ".0");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 2;
//     shape[1] = 2;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     tsr1[loc] = 0.3;  // loc is {0, 0}
//     loc[1] = 1;
//     tsr1[loc] = 0.7;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     tsr1[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
//     loc[1] = 1;
//     tsr1[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.7 * p_fail);  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.7 * p_pop);  // loc is {1, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_one_dye_first_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
    int edmans = 0;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * p_fail);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == (0.3 + 0.7) * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_one_dye_first_edman_test,
//                      *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 1;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, "0");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 2;
//     shape[1] = 2;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     tsr1[loc] = 0.3;  // loc is {0, 0}
//     loc[1] = 1;
//     tsr1[loc] = 0.7;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     tsr1[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
//     loc[1] = 1;
//     tsr1[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.7 * p_fail);  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == (0.3 + 0.7) * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {1, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_two_dyes_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 2;
    int num_channels = 1;
    DyeSeq ds(num_channels, "00");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.3;  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 0.4;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 0.5;  // loc is {1, 1}
    loc[1] = 2;
    tsr[loc] = 0.0;  // loc is {1, 2} -- one edman incompatible with 2 dyes.
    loc[0] = 2;
    loc[1] = 0;
    tsr[loc] = -1000.0;  // loc is {2, 0} -- this value should be ignored.
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {2, 1} -- this value should be ignored.
    loc[1] = 2;
    tsr[loc] = -1000.0;  // loc is {2, 2} -- this value should be ignored.
    int edmans = 1;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.1 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * p_fail);  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 0.3 * p_fail);  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    // loc is {1, 0}
    BOOST_TEST(tsr[loc] == (0.1 + 0.2 / 2.0) * p_pop + 0.4 * p_fail);
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(tsr[loc] == (0.2 / 2.0 + 0.3) * p_pop + 0.5 * p_fail);
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {1, 2}
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == (0.4 + 0.5) * p_pop);  // loc is {2, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {2, 1}
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {2, 2}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_two_dyes_second_edman_test,
//                      *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 2;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, "00");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 3;
//     shape[1] = 3;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     tsr1[loc] = 0.1;  // loc is {0, 0}
//     loc[1] = 1;
//     tsr1[loc] = 0.2;  // loc is {0, 1}
//     loc[1] = 2;
//     tsr1[loc] = 0.3;  // loc is {0, 2}
//     loc[0] = 1;
//     loc[1] = 0;
//     tsr1[loc] = 0.4;  // loc is {1, 0}
//     loc[1] = 1;
//     tsr1[loc] = 0.5;  // loc is {1, 1}
//     loc[1] = 2;
//     tsr1[loc] = 0.0;  // loc is {1, 2} -- one edman incompatible with 2 dyes.
//     loc[0] = 2;
//     loc[1] = 0;
//     tsr1[loc] = -1000.0;  // loc is {2, 0} -- this value should be ignored.
//     loc[1] = 1;
//     tsr1[loc] = -1000.0;  // loc is {2, 1} -- this value should be ignored.
//     loc[1] = 2;
//     tsr1[loc] = -1000.0;  // loc is {2, 2} -- this value should be ignored.
//     int edmans = 1;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 2);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.1 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.2 * p_fail);  // loc is {0, 1}
//     loc[1] = 2;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_fail);  // loc is {0, 2}
//     loc[0] = 1;
//     loc[1] = 0;
//     // loc is {1, 0}
//     BOOST_TEST(tsr2[loc] == (0.1 + 0.2 / 2.0) * p_pop + 0.4 * p_fail);
//     loc[1] = 1;
//     // loc is {1, 1}
//     BOOST_TEST(tsr2[loc] == (0.2 / 2.0 + 0.3) * p_pop + 0.5 * p_fail);
//     loc[1] = 2;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {1, 2}
//     loc[0] = 2;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == (0.4 + 0.5) * p_pop);  // loc is {2, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {2, 1}
//     loc[1] = 2;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {2, 2}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_three_dyes_first_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 3;
    int num_channels = 1;
    DyeSeq ds(num_channels, "000");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 4;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.3;  // loc is {0, 2}
    loc[1] = 3;
    tsr[loc] = 0.4;  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
    loc[1] = 2;
    tsr[loc] = -1000.0;  // loc is {1, 2} -- this value should be ignored.
    loc[1] = 3;
    tsr[loc] = -1000.0;  // loc is {1, 3} -- this value should be ignored.
    int edmans = 0;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.1 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * p_fail);  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 0.3 * p_fail);  // loc is {0, 2}
    loc[1] = 3;
    BOOST_TEST(tsr[loc] == 0.4 * p_fail);  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == (0.1 + 0.2 / 3.0) * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(tsr[loc] == (0.2 * 2.0 / 3.0 + 0.3 * 2.0 / 3.0) * p_pop);
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == (0.3 / 3.0 + 0.4) * p_pop);  // loc is {1, 2}
    loc[1] = 3;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {1, 3}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_three_dyes_first_edman_test,
//                      *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 3;
//     int num_channels = 1;
//     DyeSeq ds(num_channels, "000");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 2;
//     shape[1] = 4;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     tsr1[loc] = 0.1;  // loc is {0, 0}
//     loc[1] = 1;
//     tsr1[loc] = 0.2;  // loc is {0, 1}
//     loc[1] = 2;
//     tsr1[loc] = 0.3;  // loc is {0, 2}
//     loc[1] = 3;
//     tsr1[loc] = 0.4;  // loc is {0, 3}
//     loc[0] = 1;
//     loc[1] = 0;
//     tsr1[loc] = -1000.0;  // loc is {1, 0} -- this value should be ignored.
//     loc[1] = 1;
//     tsr1[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
//     loc[1] = 2;
//     tsr1[loc] = -1000.0;  // loc is {1, 2} -- this value should be ignored.
//     loc[1] = 3;
//     tsr1[loc] = -1000.0;  // loc is {1, 3} -- this value should be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == 0.1 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(tsr2[loc] == 0.2 * p_fail);  // loc is {0, 1}
//     loc[1] = 2;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_fail);  // loc is {0, 2}
//     loc[1] = 3;
//     BOOST_TEST(tsr2[loc] == 0.4 * p_fail);  // loc is {0, 3}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(tsr2[loc] == (0.1 + 0.2 / 3.0) * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     // loc is {1, 1}
//     BOOST_TEST(tsr2[loc] == (0.2 * 2.0 / 3.0 + 0.3 * 2.0 / 3.0) * p_pop);
//     loc[1] = 2;
//     BOOST_TEST(tsr2[loc] == (0.3 / 3.0 + 0.4) * p_pop);  // loc is {1, 2}
//     loc[1] = 3;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {1, 3}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_two_dye_colors_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 2;
    int num_channels = 2;
    DyeSeq ds(num_channels, "01");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.5;  // loc is {1, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.6;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 0.0;  // loc is {1, 1, 0} -- impossible state.
    loc[2] = 1;
    tsr[loc] = 0.0;  // loc is {1, 1, 1} -- impossible state.
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {2, 0, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {2, 0, 1} -- this value should be ignored.
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {2, 1, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {2, 1, 1} -- this value should be ignored.
    int edmans = 1;
    et.forward(&edmans, &tsr);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {1, 0, 0}
    BOOST_TEST(tsr[loc] == (0.1 + 0.3) * p_pop + 0.5 * p_fail);
    loc[2] = 1;
    // loc is {1, 0, 1}
    BOOST_TEST(tsr[loc] == (0.2 + 0.4) * p_pop + 0.6 * p_fail);
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {1, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {1, 1, 1}
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == (0.5 + 0.6) * p_pop);  // loc is {2, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {2, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {2, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.0);  // loc is {2, 1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_two_dye_colors_second_edman_test,
//                      *tolerance(TOL)) {
//     double p_fail = 0.05;
//     double p_pop = 0.95;
//     int num_timesteps = 2;
//     int num_channels = 2;
//     DyeSeq ds(num_channels, "01");
//     DyeTrack dt(num_timesteps, num_channels, ds);
//     EdmanTransition et(p_fail, ds, dt);
//     int order = 3;
//     int* shape = new int[order];
//     shape[0] = 3;
//     shape[1] = 2;
//     shape[2] = 2;
//     Tensor tsr1(order, shape);
//     Tensor tsr2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     tsr1[loc] = 0.1;  // loc is {0, 0, 0}
//     loc[2] = 1;
//     tsr1[loc] = 0.2;  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     tsr1[loc] = 0.3;  // loc is {0, 1, 0}
//     loc[2] = 1;
//     tsr1[loc] = 0.4;  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;
//     tsr1[loc] = 0.5;  // loc is {1, 0, 0}
//     loc[2] = 1;
//     tsr1[loc] = 0.6;  // loc is {1, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     tsr1[loc] = 0.0;  // loc is {1, 1, 0} -- impossible state.
//     loc[2] = 1;
//     tsr1[loc] = 0.0;  // loc is {1, 1, 1} -- impossible state.
//     loc[0] = 2;
//     loc[1] = 0;
//     loc[2] = 0;
//     tsr1[loc] = -1000.0;  // loc is {2, 0, 0} -- this value should be ignored.
//     loc[2] = 1;
//     tsr1[loc] = -1000.0;  // loc is {2, 0, 1} -- this value should be ignored.
//     loc[1] = 1;
//     loc[2] = 0;
//     tsr1[loc] = -1000.0;  // loc is {2, 1, 0} -- this value should be ignored.
//     loc[2] = 1;
//     tsr1[loc] = -1000.0;  // loc is {2, 1, 1} -- this value should be ignored.
//     int edmans = 1;
//     et.forward(tsr1, &edmans, &tsr2);
//     BOOST_TEST(edmans == 2);
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;
//     // loc is {1, 0, 0}
//     BOOST_TEST(tsr2[loc] == (0.1 + 0.3) * p_pop + 0.5 * p_fail);
//     loc[2] = 1;
//     // loc is {1, 0, 1}
//     BOOST_TEST(tsr2[loc] == (0.2 + 0.4) * p_pop + 0.6 * p_fail);
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {1, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {1, 1, 1}
//     loc[0] = 2;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == (0.5 + 0.6) * p_pop);  // loc is {2, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {2, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {2, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(tsr2[loc] == 0.0);  // loc is {2, 1, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(backward_in_place_trivial_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[0] = 1;
    tsr[loc] = 0.7;  // loc is {1, 0}
    int edmans = 1;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 0.7);  // loc is {0, 0}
    // We ignore the value at {1, 0}, it doesn't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_trivial_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 1;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr1[loc] = 0.3;  // loc is {0, 0}
    loc[0] = 1;
    tsr1[loc] = 0.7;  // loc is {1, 0}
    int edmans = 1;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 0.7);  // loc is {0, 0}
    // We ignore the value at {1, 0}, it doesn't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_basic_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_basic_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr1[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr1[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr1[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    tsr1[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_more_edmans_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 3;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 4;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.2;  // loc is {0, 0}
    loc[0] = 1;
    tsr[loc] = 0.3;  // loc is {1, 0}
    loc[0] = 2;
    tsr[loc] = 0.5;  // loc is {2, 0}
    loc[0] = 3;
    tsr[loc] = 0.7;  // loc is {3, 0}
    int edmans = 3;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.2 + p_pop * 0.3);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {1, 0}
    loc[0] = 2;
    BOOST_TEST(tsr[loc] == p_fail * 0.5 + p_pop * 0.7);  // loc is {2, 0}
    // We ignore the value at {3, 0}, it doesn't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_more_edmans_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 3;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 4;
    shape[1] = 1;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr1[loc] = 0.2;  // loc is {0, 0}
    loc[0] = 1;
    tsr1[loc] = 0.3;  // loc is {1, 0}
    loc[0] = 2;
    tsr1[loc] = 0.5;  // loc is {2, 0}
    loc[0] = 3;
    tsr1[loc] = 0.7;  // loc is {3, 0}
    int edmans = 3;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.2 + p_pop * 0.3);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {1, 0}
    loc[0] = 2;
    BOOST_TEST(tsr2[loc] == p_fail * 0.5 + p_pop * 0.7);  // loc is {2, 0}
    // We ignore the value at {3, 0}, it doesn't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 1.11;  // loc is {1, 0, 0}
    loc[2] = 1;
    tsr[loc] = 1.22;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 1.33;  // loc is {1, 1, 0}
    loc[2] = 1;
    tsr[loc] = 1.44;  // loc is {1, 1, 1}
    int edmans = 1;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.2 + p_pop * 1.22);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.4 + p_pop * 1.44);  // loc is {0, 1, 1}
    // We ignore the values at {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, and {1, 1, 1}
    // they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr1[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr1[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr1[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr1[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    tsr1[loc] = 1.11;  // loc is {1, 0, 0}
    loc[2] = 1;
    tsr1[loc] = 1.22;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr1[loc] = 1.33;  // loc is {1, 1, 0}
    loc[2] = 1;
    tsr1[loc] = 1.44;  // loc is {1, 1, 1}
    int edmans = 1;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.2 + p_pop * 1.22);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.4 + p_pop * 1.44);  // loc is {0, 1, 1}
    // We ignore the values at {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, and {1, 1, 1}
    // they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_irrelevant_dye_seq_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, ".0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_irrelevant_dye_seq_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, ".0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr1[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr1[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr1[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    tsr1[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_one_dye_first_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
    int edmans = 1;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.7 + p_pop * 1.33);  // loc is {0, 1}
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_one_dye_first_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr1[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr1[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr1[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    tsr1[loc] = -1000.0;  // loc is {1, 1} -- this value should be ignored.
    int edmans = 1;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.7 + p_pop * 1.33);  // loc is {0, 1}
    // We ignore the values at {1, 0} and {1, 1}, they don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_two_dyes_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 2;
    int num_channels = 1;
    DyeSeq ds(num_channels, "00");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.3;  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 0.4;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 0.5;  // loc is {1, 1}
    loc[1] = 2;
    tsr[loc] = -1000.0;  // loc is {1, 2} -- this value should be ignored.
    loc[0] = 2;
    loc[1] = 0;
    tsr[loc] = 0.6;  // loc is {2, 0}
    loc[1] = 1;
    tsr[loc] = -1000.0;  // loc is {2, 1} -- this value should be ignored.
    loc[1] = 2;
    tsr[loc] = -1000.0;  // loc is {2, 2} -- this value should be ignored.
    int edmans = 2;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.1 + p_pop * 0.4);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(tsr[loc] == p_fail * 0.2 + p_pop * (0.5 / 2.0 + 0.4 / 2.0));
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    // loc is {1, 0}
    BOOST_TEST(tsr[loc] == p_fail * 0.4 + p_pop * 0.6);
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(tsr[loc] == p_fail * 0.5 + p_pop * 0.6);
    // We ignore the values at {1, 2}, {2, 0}, {2, 1}, and {2, 2}.
    // They don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_two_dyes_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 2;
    int num_channels = 1;
    DyeSeq ds(num_channels, "00");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 3;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr1[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    tsr1[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    tsr1[loc] = 0.3;  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    tsr1[loc] = 0.4;  // loc is {1, 0}
    loc[1] = 1;
    tsr1[loc] = 0.5;  // loc is {1, 1}
    loc[1] = 2;
    tsr1[loc] = -1000.0;  // loc is {1, 2} -- this value should be ignored.
    loc[0] = 2;
    loc[1] = 0;
    tsr1[loc] = 0.6;  // loc is {2, 0}
    loc[1] = 1;
    tsr1[loc] = -1000.0;  // loc is {2, 1} -- this value should be ignored.
    loc[1] = 2;
    tsr1[loc] = -1000.0;  // loc is {2, 2} -- this value should be ignored.
    int edmans = 2;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.1 + p_pop * 0.4);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(tsr2[loc] == p_fail * 0.2 + p_pop * (0.5 / 2.0 + 0.4 / 2.0));
    loc[1] = 2;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    // loc is {1, 0}
    BOOST_TEST(tsr2[loc] == p_fail * 0.4 + p_pop * 0.6);
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(tsr2[loc] == p_fail * 0.5 + p_pop * 0.6);
    // We ignore the values at {1, 2}, {2, 0}, {2, 1}, and {2, 2}.
    // They don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_three_dyes_first_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 3;
    int num_channels = 1;
    DyeSeq ds(num_channels, "000");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 4;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.3;  // loc is {0, 2}
    loc[1] = 3;
    tsr[loc] = 0.4;  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 1.11;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 1.22;  // loc is {1, 1}
    loc[1] = 2;
    tsr[loc] = 1.33;  // loc is {1, 2}
    loc[1] = 3;
    tsr[loc] = -1000.0;  // loc is {1, 3} -- this value should be ignored.
    int edmans = 1;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(tsr[loc]
               == p_fail * 0.2 + p_pop * (1.11 * 1.0 / 3.0 + 1.22 * 2.0 / 3.0));
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(tsr[loc]
               == p_fail * 0.3 + p_pop * (1.22 * 2.0 / 3.0 + 1.33 * 1.0 / 3.0));
    loc[1] = 3;
    BOOST_TEST(tsr[loc] == p_fail * 0.4 + p_pop * 1.33);  // loc is {0, 3}
    // We ignore the values at {1, 0}, {1, 1}, {1, 2}, and {1, 3}.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_three_dyes_first_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 3;
    int num_channels = 1;
    DyeSeq ds(num_channels, "000");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 4;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr1[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    tsr1[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    tsr1[loc] = 0.3;  // loc is {0, 2}
    loc[1] = 3;
    tsr1[loc] = 0.4;  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    tsr1[loc] = 1.11;  // loc is {1, 0}
    loc[1] = 1;
    tsr1[loc] = 1.22;  // loc is {1, 1}
    loc[1] = 2;
    tsr1[loc] = 1.33;  // loc is {1, 2}
    loc[1] = 3;
    tsr1[loc] = -1000.0;  // loc is {1, 3} -- this value should be ignored.
    int edmans = 1;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(tsr2[loc]
               == p_fail * 0.2 + p_pop * (1.11 * 1.0 / 3.0 + 1.22 * 2.0 / 3.0));
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(tsr2[loc]
               == p_fail * 0.3 + p_pop * (1.22 * 2.0 / 3.0 + 1.33 * 1.0 / 3.0));
    loc[1] = 3;
    BOOST_TEST(tsr2[loc] == p_fail * 0.4 + p_pop * 1.33);  // loc is {0, 3}
    // We ignore the values at {1, 0}, {1, 1}, {1, 2}, and {1, 3}.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_two_dye_colors_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 2;
    int num_channels = 2;
    DyeSeq ds(num_channels, "01");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.5;  // loc is {1, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.6;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {1, 1, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {1, 1, 1} -- this value should be ignored.
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.7;  // loc is {2, 0, 0}
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {2, 0, 1} -- this value should be ignored.
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = -1000.0;  // loc is {2, 1, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr[loc] = -1000.0;  // loc is {2, 1, 1} -- this value should be ignored.
    int edmans = 2;
    et.backward(tsr, &edmans, &tsr);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.1 + p_pop * 0.5);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.2 + p_pop * 0.6);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == p_fail * 0.4 + p_pop * 0.6);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {1, 0, 0}
    BOOST_TEST(tsr[loc] == p_fail * 0.5 + p_pop * 0.7);
    loc[2] = 1;
    // loc is {1, 0, 1}
    BOOST_TEST(tsr[loc] == p_fail * 0.6 + p_pop * 0.7);
    // We ignore the values at {1, 1, 0}, {1, 1, 1}, {2, 0, 0}, {2, 0, 1},
    // {2, 1, 0}, and {2, 1, 1}. They don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_two_dye_colors_second_edman_test,
                     *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 2;
    int num_channels = 2;
    DyeSeq ds(num_channels, "01");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr1[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr1[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr1[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr1[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    tsr1[loc] = 0.5;  // loc is {1, 0, 0}
    loc[2] = 1;
    tsr1[loc] = 0.6;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr1[loc] = -1000.0;  // loc is {1, 1, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr1[loc] = -1000.0;  // loc is {1, 1, 1} -- this value should be ignored.
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    tsr1[loc] = 0.7;  // loc is {2, 0, 0}
    loc[2] = 1;
    tsr1[loc] = -1000.0;  // loc is {2, 0, 1} -- this value should be ignored.
    loc[1] = 1;
    loc[2] = 0;
    tsr1[loc] = -1000.0;  // loc is {2, 1, 0} -- this value should be ignored.
    loc[2] = 1;
    tsr1[loc] = -1000.0;  // loc is {2, 1, 1} -- this value should be ignored.
    int edmans = 2;
    et.backward(tsr1, &edmans, &tsr2);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.1 + p_pop * 0.5);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.2 + p_pop * 0.6);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == p_fail * 0.4 + p_pop * 0.6);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {1, 0, 0}
    BOOST_TEST(tsr2[loc] == p_fail * 0.5 + p_pop * 0.7);
    loc[2] = 1;
    // loc is {1, 0, 1}
    BOOST_TEST(tsr2[loc] == p_fail * 0.6 + p_pop * 0.7);
    // We ignore the values at {1, 1, 0}, {1, 1, 1}, {2, 0, 0}, {2, 0, 1},
    // {2, 1, 0}, and {2, 1, 1}. They don't matter.
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor ftsr(order, shape);
    Tensor btsr(order, shape);
    Tensor nbtsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    ftsr[loc] = 0.61;  // loc is {0, 0}
    btsr[loc] = 0.51;
    nbtsr[loc] = 0.41;
    loc[1] = 1;
    ftsr[loc] = 0.91;  // loc is {0, 1}
    btsr[loc] = 0.81;
    nbtsr[loc] = 0.71;
    loc[0] = 1;
    loc[1] = 0;
    ftsr[loc] = 0.62;  // loc is {1, 0}
    btsr[loc] = 0.52;
    nbtsr[loc] = 0.42;
    loc[1] = 1;
    ftsr[loc] = 0.92;  // loc is {1, 1}
    btsr[loc] = 0.82;
    nbtsr[loc] = 0.72;
    delete[] loc;
    int edmans = 0;
    double probability = 1.0;
    ErrorModelFitter emf(DistributionType::LOGNORMAL);
    et.improve_fit(ftsr, btsr, nbtsr, edmans, probability, &emf);
    BOOST_TEST(emf.p_edman_failure_fit.get()
               == (0.91 * p_fail * 0.71) / (0.91 * 0.81));
}

BOOST_AUTO_TEST_CASE(improve_fit_twice_test, *tolerance(TOL)) {
    double p_fail = 0.05;
    double p_pop = 0.95;
    int num_timesteps = 1;
    int num_channels = 1;
    DyeSeq ds(num_channels, "");
    DyeTrack dt(num_timesteps, num_channels, ds);
    EdmanTransition et(p_fail, ds, dt);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    Tensor ftsr1(order, shape);
    Tensor btsr1(order, shape);
    Tensor nbtsr1(order, shape);
    Tensor ftsr2(order, shape);
    Tensor btsr2(order, shape);
    Tensor nbtsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    ftsr1[loc] = 0.31;  // loc is {0, 0}
    btsr1[loc] = 0.331;
    nbtsr1[loc] = 0.21;
    ftsr2[loc] = 0.221;
    btsr2[loc] = 0.11;
    nbtsr2[loc] = 0.111;
    loc[1] = 1;
    ftsr1[loc] = 0.91;  // loc is {0, 1}
    btsr1[loc] = 0.81;
    nbtsr1[loc] = 0.71;
    ftsr2[loc] = 0.61;
    btsr2[loc] = 0.51;
    nbtsr2[loc] = 0.41;
    loc[0] = 1;
    loc[1] = 0;
    ftsr1[loc] = 0.32;  // loc is {1, 0}
    btsr1[loc] = 0.332;
    nbtsr1[loc] = 0.22;
    ftsr2[loc] = 0.222;
    btsr2[loc] = 0.12;
    nbtsr2[loc] = 0.112;
    loc[1] = 1;
    ftsr1[loc] = 0.92;  // loc is {1, 1}
    btsr1[loc] = 0.82;
    nbtsr1[loc] = 0.72;
    ftsr2[loc] = 0.92;
    btsr2[loc] = 0.82;
    nbtsr2[loc] = 0.72;
    delete[] loc;
    int edmans = 0;
    double prob1 = 0.12345;
    double prob2 = 0.98765;
    ErrorModelFitter emf(DistributionType::LOGNORMAL);
    et.improve_fit(ftsr1, btsr1, nbtsr1, edmans, prob1, &emf);
    et.improve_fit(ftsr2, btsr2, nbtsr2, edmans, prob2, &emf);
    BOOST_TEST(emf.p_edman_failure_fit.get()
               == (0.91 * p_fail * 0.71 / prob1 + 0.61 * p_fail * 0.41 / prob2)
                          / (0.91 * 0.81 / prob1 + 0.61 * 0.51 / prob2));
}

BOOST_AUTO_TEST_SUITE_END()  // edman_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
