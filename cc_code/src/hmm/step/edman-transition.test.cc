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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 1.0;  // loc is {0, 0}
    loc[0] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
    int edmans = 0;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 1.0 * p_fail);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(psv.tensor[loc] == 1.0 * p_pop);  // loc is {1, 0}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 1.0;  // loc is {0, 0}
//     loc[0] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 1.0 * p_fail);  // loc is {0, 0}
//     loc[0] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 1.0 * p_pop);  // loc is {1, 0}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
    loc[1] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
    int edmans = 0;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p_fail);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p_pop);  // loc is {1, 1}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.7;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
//     loc[1] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p_fail);  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p_pop);  // loc is {1, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[0] = 1;
    psv.tensor[loc] = 0.3;  // loc is {1, 0}
    loc[0] = 2;
    psv.tensor[loc] = 0.5;  // loc is {2, 0}
    loc[0] = 3;
    psv.tensor[loc] = -1000.0;  // loc is {3, 0} -- to be ignored.
    int edmans = 2;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 3);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.2 * p_fail);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * p_pop + 0.3 * p_fail);  // loc is {1, 0}
    loc[0] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_pop + 0.5 * p_fail);  // loc is {2, 0}
    loc[0] = 3;
    BOOST_TEST(psv.tensor[loc] == 0.5 * p_pop);  // loc is {3, 0}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.2;  // loc is {0, 0}
//     loc[0] = 1;
//     psv1.tensor[loc] = 0.3;  // loc is {1, 0}
//     loc[0] = 2;
//     psv1.tensor[loc] = 0.5;  // loc is {2, 0}
//     loc[0] = 3;
//     psv1.tensor[loc] = -1000.0;  // loc is {3, 0} -- to be ignored.
//     int edmans = 2;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 3);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 * p_fail);  // loc is {0, 0}
//     loc[0] = 1;
//     BOOST_TEST(psv2.tensor[loc]
//                == 0.2 * p_pop + 0.3 * p_fail);  // loc is {1, 0}
//     loc[0] = 2;
//     BOOST_TEST(psv2.tensor[loc]
//                == 0.3 * p_pop + 0.5 * p_fail);  // loc is {2, 0}
//     loc[0] = 3;
//     BOOST_TEST(psv2.tensor[loc]
//                == 0.5 * p_pop);  // loc is {3, 0}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {1, 0, 0} -- to be ignored.
    loc[2] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 0, 1} -- to be ignored.
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1, 0} -- to be ignored.
    loc[2] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1, 1} -- to be ignored.
    int edmans = 0;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1 * p_pop);  // loc is {1, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * p_pop);  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_pop);  // loc is {1, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.4 * p_pop);  // loc is {1, 1, 1}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     psv1.tensor[loc] = 0.1;  // loc is {0, 0, 0}
//     loc[2] = 1;
//     psv1.tensor[loc] = 0.2;  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 1, 0}
//     loc[2] = 1;
//     psv1.tensor[loc] = 0.4;  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 0, 0} -- ignore this
//     loc[2] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 0, 1} -- ignore this
//     loc[1] = 1;
//     loc[2] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 1, 0} -- to be ignored.
//     loc[2] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 1, 1} -- to be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.1 * p_pop);  // loc is {1, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 * p_pop);  // loc is {1, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_pop);  // loc is {1, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.4 * p_pop);  // loc is {1, 1, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
    loc[1] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
    int edmans = 0;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p_fail);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p_pop);  // loc is {1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_irrelevant_dye_seq_test,
// *tolerance(TOL)) {
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.7;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
//     loc[1] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p_fail);  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p_pop);  // loc is {1, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
    loc[1] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
    int edmans = 0;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p_fail);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == (0.3 + 0.7) * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {1, 1}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.7;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
//     loc[1] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p_fail);  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == (0.3 + 0.7) * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {1, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.3;  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 0.4;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.5;  // loc is {1, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.0;  // loc is {1, 2} -- 1 edman can't have 2 dyes.
    loc[0] = 2;
    loc[1] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {2, 0} -- to be ignored.
    loc[1] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {2, 1} -- to be ignored.
    loc[1] = 2;
    psv.tensor[loc] = -1000.0;  // loc is {2, 2} -- to be ignored.
    int edmans = 1;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * p_fail);  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_fail);  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    // loc is {1, 0}
    BOOST_TEST(psv.tensor[loc] == (0.1 + 0.2 / 2.0) * p_pop + 0.4 * p_fail);
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(psv.tensor[loc] == (0.2 / 2.0 + 0.3) * p_pop + 0.5 * p_fail);
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {1, 2}
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == (0.4 + 0.5) * p_pop);  // loc is {2, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {2, 1}
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {2, 2}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.1;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.2;  // loc is {0, 1}
//     loc[1] = 2;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 2}
//     loc[0] = 1;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.4;  // loc is {1, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.5;  // loc is {1, 1}
//     loc[1] = 2;
//     psv1.tensor[loc] = 0.0;  // loc is {1, 2} -- 1 edman can't have 2 dyes.
//     loc[0] = 2;
//     loc[1] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {2, 0} -- to be ignored.
//     loc[1] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {2, 1} -- to be ignored.
//     loc[1] = 2;
//     psv1.tensor[loc] = -1000.0;  // loc is {2, 2} -- to be ignored.
//     int edmans = 1;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 2);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 * p_fail);  // loc is {0, 1}
//     loc[1] = 2;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_fail);  // loc is {0, 2}
//     loc[0] = 1;
//     loc[1] = 0;
//     // loc is {1, 0}
//     BOOST_TEST(psv2.tensor[loc] == (0.1 + 0.2 / 2.0) * p_pop + 0.4 * p_fail);
//     loc[1] = 1;
//     // loc is {1, 1}
//     BOOST_TEST(psv2.tensor[loc] == (0.2 / 2.0 + 0.3) * p_pop + 0.5 * p_fail);
//     loc[1] = 2;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {1, 2}
//     loc[0] = 2;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == (0.4 + 0.5) * p_pop);  // loc is {2, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {2, 1}
//     loc[1] = 2;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {2, 2}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.3;  // loc is {0, 2}
    loc[1] = 3;
    psv.tensor[loc] = 0.4;  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
    loc[1] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
    loc[1] = 2;
    psv.tensor[loc] = -1000.0;  // loc is {1, 2} -- to be ignored.
    loc[1] = 3;
    psv.tensor[loc] = -1000.0;  // loc is {1, 3} -- to be ignored.
    int edmans = 0;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * p_fail);  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_fail);  // loc is {0, 2}
    loc[1] = 3;
    BOOST_TEST(psv.tensor[loc] == 0.4 * p_fail);  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == (0.1 + 0.2 / 3.0) * p_pop);  // loc is {1, 0}
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(psv.tensor[loc] == (0.2 * 2.0 / 3.0 + 0.3 * 2.0 / 3.0) * p_pop);
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == (0.3 / 3.0 + 0.4) * p_pop);  // loc is {1, 2}
    loc[1] = 3;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {1, 3}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.1;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.2;  // loc is {0, 1}
//     loc[1] = 2;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 2}
//     loc[1] = 3;
//     psv1.tensor[loc] = 0.4;  // loc is {0, 3}
//     loc[0] = 1;
//     loc[1] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 0} -- to be ignored.
//     loc[1] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
//     loc[1] = 2;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 2} -- to be ignored.
//     loc[1] = 3;
//     psv1.tensor[loc] = -1000.0;  // loc is {1, 3} -- to be ignored.
//     int edmans = 0;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 1);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 * p_fail);  // loc is {0, 1}
//     loc[1] = 2;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_fail);  // loc is {0, 2}
//     loc[1] = 3;
//     BOOST_TEST(psv2.tensor[loc] == 0.4 * p_fail);  // loc is {0, 3}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc]
//                == (0.1 + 0.2 / 3.0) * p_pop);  // loc is {1, 0}
//     loc[1] = 1;
//     // loc is {1, 1}
//     BOOST_TEST(psv2.tensor[loc]
//                == (0.2 * 2.0 / 3.0 + 0.3 * 2.0 / 3.0) * p_pop);
//     loc[1] = 2;
//     BOOST_TEST(psv2.tensor[loc]
//                == (0.3 / 3.0 + 0.4) * p_pop);  // loc is {1, 2}
//     loc[1] = 3;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {1, 3}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 0.5;  // loc is {1, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.6;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = 0.0;  // loc is {1, 1, 0} -- impossible state.
    loc[2] = 1;
    psv.tensor[loc] = 0.0;  // loc is {1, 1, 1} -- impossible state.
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {2, 0, 0} -- to be ignored.
    loc[2] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {2, 0, 1} -- to be ignored.
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {2, 1, 0} -- to be ignored.
    loc[2] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {2, 1, 1} -- to be ignored.
    int edmans = 1;
    et.forward(&edmans, &psv);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {1, 0, 0}
    BOOST_TEST(psv.tensor[loc] == (0.1 + 0.3) * p_pop + 0.5 * p_fail);
    loc[2] = 1;
    // loc is {1, 0, 1}
    BOOST_TEST(psv.tensor[loc] == (0.2 + 0.4) * p_pop + 0.6 * p_fail);
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {1, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {1, 1, 1}
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == (0.5 + 0.6) * p_pop);  // loc is {2, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {2, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {2, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.0);  // loc is {2, 1, 1}
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
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     psv1.tensor[loc] = 0.1;  // loc is {0, 0, 0}
//     loc[2] = 1;
//     psv1.tensor[loc] = 0.2;  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 1, 0}
//     loc[2] = 1;
//     psv1.tensor[loc] = 0.4;  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;
//     psv1.tensor[loc] = 0.5;  // loc is {1, 0, 0}
//     loc[2] = 1;
//     psv1.tensor[loc] = 0.6;  // loc is {1, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     psv1.tensor[loc] = 0.0;  // loc is {1, 1, 0} -- impossible state.
//     loc[2] = 1;
//     psv1.tensor[loc] = 0.0;  // loc is {1, 1, 1} -- impossible state.
//     loc[0] = 2;
//     loc[1] = 0;
//     loc[2] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {2, 0, 0} -- to be ignored
//     loc[2] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {2, 0, 1} -- to be ignored
//     loc[1] = 1;
//     loc[2] = 0;
//     psv1.tensor[loc] = -1000.0;  // loc is {2, 1, 0} -- to be ignored
//     loc[2] = 1;
//     psv1.tensor[loc] = -1000.0;  // loc is {2, 1, 1} -- to be ignored.
//     int edmans = 1;
//     et.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(edmans == 2);
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.1 * p_fail);  // loc is {0, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 * p_fail);  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p_fail);  // loc is {0, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.4 * p_fail);  // loc is {0, 1, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     loc[2] = 0;  // loc is {1, 0, 0}
//     BOOST_TEST(psv2.tensor[loc] == (0.1 + 0.3) * p_pop + 0.5 * p_fail);
//     loc[2] = 1;
//     // loc is {1, 0, 1}
//     BOOST_TEST(psv2.tensor[loc] == (0.2 + 0.4) * p_pop + 0.6 * p_fail);
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {1, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {1, 1, 1}
//     loc[0] = 2;
//     loc[1] = 0;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == (0.5 + 0.6) * p_pop);  // loc is {2, 0, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {2, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {2, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.0);  // loc is {2, 1, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[0] = 1;
    psv.tensor[loc] = 0.7;  // loc is {1, 0}
    int edmans = 1;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.3 + p_pop * 0.7);  // loc is {0, 0}
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[0] = 1;
    psv1.tensor[loc] = 0.7;  // loc is {1, 0}
    int edmans = 1;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 0.7);  // loc is {0, 0}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv1.tensor[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[0] = 1;
    psv.tensor[loc] = 0.3;  // loc is {1, 0}
    loc[0] = 2;
    psv.tensor[loc] = 0.5;  // loc is {2, 0}
    loc[0] = 3;
    psv.tensor[loc] = 0.7;  // loc is {3, 0}
    int edmans = 3;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.2 + p_pop * 0.3);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {1, 0}
    loc[0] = 2;
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.5 + p_pop * 0.7);  // loc is {2, 0}
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[0] = 1;
    psv1.tensor[loc] = 0.3;  // loc is {1, 0}
    loc[0] = 2;
    psv1.tensor[loc] = 0.5;  // loc is {2, 0}
    loc[0] = 3;
    psv1.tensor[loc] = 0.7;  // loc is {3, 0}
    int edmans = 3;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.2 + p_pop * 0.3);  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 0.5);  // loc is {1, 0}
    loc[0] = 2;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.5 + p_pop * 0.7);  // loc is {2, 0}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 1.11;  // loc is {1, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = 1.22;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = 1.33;  // loc is {1, 1, 0}
    loc[2] = 1;
    psv.tensor[loc] = 1.44;  // loc is {1, 1, 1}
    int edmans = 1;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.2 + p_pop * 1.22);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.4 + p_pop * 1.44);  // loc is {0, 1, 1}
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    psv1.tensor[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    psv1.tensor[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    psv1.tensor[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    psv1.tensor[loc] = 1.11;  // loc is {1, 0, 0}
    loc[2] = 1;
    psv1.tensor[loc] = 1.22;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv1.tensor[loc] = 1.33;  // loc is {1, 1, 0}
    loc[2] = 1;
    psv1.tensor[loc] = 1.44;  // loc is {1, 1, 1}
    int edmans = 1;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.2 + p_pop * 1.22);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.4 + p_pop * 1.44);  // loc is {0, 1, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv1.tensor[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 1.77;  // loc is {1, 1}
    int edmans = 1;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.7 + p_pop * 1.77);  // loc is {0, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
    int edmans = 1;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.7 + p_pop * 1.33);  // loc is {0, 1}
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.7;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv1.tensor[loc] = 1.33;  // loc is {1, 0}
    loc[1] = 1;
    psv1.tensor[loc] = -1000.0;  // loc is {1, 1} -- to be ignored.
    int edmans = 1;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 1.33);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.7 + p_pop * 1.33);  // loc is {0, 1}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.3;  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 0.4;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.5;  // loc is {1, 1}
    loc[1] = 2;
    psv.tensor[loc] = -1000.0;  // loc is {1, 2} -- to be ignored.
    loc[0] = 2;
    loc[1] = 0;
    psv.tensor[loc] = 0.6;  // loc is {2, 0}
    loc[1] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {2, 1} -- to be ignored.
    loc[1] = 2;
    psv.tensor[loc] = -1000.0;  // loc is {2, 2} -- to be ignored.
    int edmans = 2;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.1 + p_pop * 0.4);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.2 + p_pop * (0.5 / 2.0 + 0.4 / 2.0));
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    // loc is {1, 0}
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.4 + p_pop * 0.6);
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.5 + p_pop * 0.6);
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    psv1.tensor[loc] = 0.3;  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    psv1.tensor[loc] = 0.4;  // loc is {1, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.5;  // loc is {1, 1}
    loc[1] = 2;
    psv1.tensor[loc] = -1000.0;  // loc is {1, 2} -- to be ignored.
    loc[0] = 2;
    loc[1] = 0;
    psv1.tensor[loc] = 0.6;  // loc is {2, 0}
    loc[1] = 1;
    psv1.tensor[loc] = -1000.0;  // loc is {2, 1} -- to be ignored.
    loc[1] = 2;
    psv1.tensor[loc] = -1000.0;  // loc is {2, 2} -- to be ignored.
    int edmans = 2;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.1 + p_pop * 0.4);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.2 + p_pop * (0.5 / 2.0 + 0.4 / 2.0));
    loc[1] = 2;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 2}
    loc[0] = 1;
    loc[1] = 0;
    // loc is {1, 0}
    BOOST_TEST(psv2.tensor[loc] == p_fail * 0.4 + p_pop * 0.6);
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(psv2.tensor[loc] == p_fail * 0.5 + p_pop * 0.6);
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.3;  // loc is {0, 2}
    loc[1] = 3;
    psv.tensor[loc] = 0.4;  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 1.11;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 1.22;  // loc is {1, 1}
    loc[1] = 2;
    psv.tensor[loc] = 1.33;  // loc is {1, 2}
    loc[1] = 3;
    psv.tensor[loc] = -1000.0;  // loc is {1, 3} -- to be ignored.
    int edmans = 1;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.2 + p_pop * (1.11 * 1.0 / 3.0 + 1.22 * 2.0 / 3.0));
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.3 + p_pop * (1.22 * 2.0 / 3.0 + 1.33 * 1.0 / 3.0));
    loc[1] = 3;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.4 + p_pop * 1.33);  // loc is {0, 3}
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[1] = 2;
    psv1.tensor[loc] = 0.3;  // loc is {0, 2}
    loc[1] = 3;
    psv1.tensor[loc] = 0.4;  // loc is {0, 3}
    loc[0] = 1;
    loc[1] = 0;
    psv1.tensor[loc] = 1.11;  // loc is {1, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 1.22;  // loc is {1, 1}
    loc[1] = 2;
    psv1.tensor[loc] = 1.33;  // loc is {1, 2}
    loc[1] = 3;
    psv1.tensor[loc] = -1000.0;  // loc is {1, 3} -- to be ignored.
    int edmans = 1;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 0);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.1 + p_pop * 1.11);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.2 + p_pop * (1.11 * 1.0 / 3.0 + 1.22 * 2.0 / 3.0));
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * (1.22 * 2.0 / 3.0 + 1.33 * 1.0 / 3.0));
    loc[1] = 3;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.4 + p_pop * 1.33);  // loc is {0, 3}
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
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 0.5;  // loc is {1, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = 0.6;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1, 0} -- to be ignored.
    loc[2] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {1, 1, 1} -- to be ignored.
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    psv.tensor[loc] = 0.7;  // loc is {2, 0, 0}
    loc[2] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {2, 0, 1} -- to be ignored.
    loc[1] = 1;
    loc[2] = 0;
    psv.tensor[loc] = -1000.0;  // loc is {2, 1, 0} -- to be ignored.
    loc[2] = 1;
    psv.tensor[loc] = -1000.0;  // loc is {2, 1, 1} -- to be ignored.
    int edmans = 2;
    et.backward(psv, &edmans, &psv);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.1 + p_pop * 0.5);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.2 + p_pop * 0.6);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc]
               == p_fail * 0.4 + p_pop * 0.6);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {1, 0, 0}
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.5 + p_pop * 0.7);
    loc[2] = 1;
    // loc is {1, 0, 1}
    BOOST_TEST(psv.tensor[loc] == p_fail * 0.6 + p_pop * 0.7);
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
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    psv1.tensor[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    psv1.tensor[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    psv1.tensor[loc] = 0.4;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    psv1.tensor[loc] = 0.5;  // loc is {1, 0, 0}
    loc[2] = 1;
    psv1.tensor[loc] = 0.6;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    psv1.tensor[loc] = -1000.0;  // loc is {1, 1, 0} -- to be ignored.
    loc[2] = 1;
    psv1.tensor[loc] = -1000.0;  // loc is {1, 1, 1} -- to be ignored.
    loc[0] = 2;
    loc[1] = 0;
    loc[2] = 0;
    psv1.tensor[loc] = 0.7;  // loc is {2, 0, 0}
    loc[2] = 1;
    psv1.tensor[loc] = -1000.0;  // loc is {2, 0, 1} -- to be ignored.
    loc[1] = 1;
    loc[2] = 0;
    psv1.tensor[loc] = -1000.0;  // loc is {2, 1, 0} -- to be ignored.
    loc[2] = 1;
    psv1.tensor[loc] = -1000.0;  // loc is {2, 1, 1} -- to be ignored.
    int edmans = 2;
    et.backward(psv1, &edmans, &psv2);
    BOOST_TEST(edmans == 1);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.1 + p_pop * 0.5);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.2 + p_pop * 0.6);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.3 + p_pop * 0.5);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv2.tensor[loc]
               == p_fail * 0.4 + p_pop * 0.6);  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {1, 0, 0}
    BOOST_TEST(psv2.tensor[loc] == p_fail * 0.5 + p_pop * 0.7);
    loc[2] = 1;
    // loc is {1, 0, 1}
    BOOST_TEST(psv2.tensor[loc] == p_fail * 0.6 + p_pop * 0.7);
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
    PeptideStateVector fpsv(order, shape);
    PeptideStateVector bpsv(order, shape);
    PeptideStateVector nbpsv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    fpsv.tensor[loc] = 0.61;  // loc is {0, 0}
    bpsv.tensor[loc] = 0.51;
    nbpsv.tensor[loc] = 0.41;
    loc[1] = 1;
    fpsv.tensor[loc] = 0.91;  // loc is {0, 1}
    bpsv.tensor[loc] = 0.81;
    nbpsv.tensor[loc] = 0.71;
    loc[0] = 1;
    loc[1] = 0;
    fpsv.tensor[loc] = 0.62;  // loc is {1, 0}
    bpsv.tensor[loc] = 0.52;
    nbpsv.tensor[loc] = 0.42;
    loc[1] = 1;
    fpsv.tensor[loc] = 0.92;  // loc is {1, 1}
    bpsv.tensor[loc] = 0.82;
    nbpsv.tensor[loc] = 0.72;
    delete[] loc;
    int edmans = 0;
    double probability = 1.0;
    SequencingModelFitter smf;
    et.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    BOOST_TEST(smf.p_edman_failure_fit.get()
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
    PeptideStateVector fpsv1(order, shape);
    PeptideStateVector bpsv1(order, shape);
    PeptideStateVector nbpsv1(order, shape);
    PeptideStateVector fpsv2(order, shape);
    PeptideStateVector bpsv2(order, shape);
    PeptideStateVector nbpsv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    fpsv1.tensor[loc] = 0.31;  // loc is {0, 0}
    bpsv1.tensor[loc] = 0.331;
    nbpsv1.tensor[loc] = 0.21;
    fpsv2.tensor[loc] = 0.221;
    bpsv2.tensor[loc] = 0.11;
    nbpsv2.tensor[loc] = 0.111;
    loc[1] = 1;
    fpsv1.tensor[loc] = 0.91;  // loc is {0, 1}
    bpsv1.tensor[loc] = 0.81;
    nbpsv1.tensor[loc] = 0.71;
    fpsv2.tensor[loc] = 0.61;
    bpsv2.tensor[loc] = 0.51;
    nbpsv2.tensor[loc] = 0.41;
    loc[0] = 1;
    loc[1] = 0;
    fpsv1.tensor[loc] = 0.32;  // loc is {1, 0}
    bpsv1.tensor[loc] = 0.332;
    nbpsv1.tensor[loc] = 0.22;
    fpsv2.tensor[loc] = 0.222;
    bpsv2.tensor[loc] = 0.12;
    nbpsv2.tensor[loc] = 0.112;
    loc[1] = 1;
    fpsv1.tensor[loc] = 0.92;  // loc is {1, 1}
    bpsv1.tensor[loc] = 0.82;
    nbpsv1.tensor[loc] = 0.72;
    fpsv2.tensor[loc] = 0.92;
    bpsv2.tensor[loc] = 0.82;
    nbpsv2.tensor[loc] = 0.72;
    delete[] loc;
    int edmans = 0;
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
