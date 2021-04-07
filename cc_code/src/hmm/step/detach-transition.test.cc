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
#include "detach-transition.h"
#include "hmm/state-vector/peptide-state-vector.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
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

BOOST_AUTO_TEST_CASE(forward_in_place_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 1.0;  // loc is {0, 0}
    int edmans = 0;
    dt.forward(&edmans, &psv);
    BOOST_TEST(psv.tensor[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_trivial_test, *tolerance(TOL)) {
//     double p_detach = 0.05;
//     DetachTransition dt(p_detach);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 1;
//     shape[1] = 1;
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 1.0;  // loc is {0, 0}
//     int edmans = 0;
//     dt.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(psv2.tensor[loc] == 1.0);  // loc is {0, 0}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    int edmans = 0;
    dt.forward(&edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 + 0.7 * p_detach);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * (1 - p_detach));  // loc is {0, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_basic_test, *tolerance(TOL)) {
//     double p_detach = 0.05;
//     DetachTransition dt(p_detach);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 1;
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
//     int edmans = 0;
//     dt.forward(tsr1, &edmans, &psv2);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 + 0.7 * p_detach);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * (1 - p_detach));  // loc is {0, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.6;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.1;  // loc is {0, 2}
    int edmans = 0;
    dt.forward(&edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc]
               == 0.3 + (0.6 + 0.1) * p_detach);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.6 * (1 - p_detach));  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.1 * (1 - p_detach));  // loc is {0, 2}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_bigger_test, *tolerance(TOL)) {
//     double p_detach = 0.05;
//     DetachTransition dt(p_detach);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 1;
//     shape[1] = 3;
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.6;  // loc is {0, 1}
//     loc[1] = 2;
//     psv1.tensor[loc] = 0.1;  // loc is {0, 2}
//     int edmans = 0;
//     dt.forward(tsr1, &edmans, &psv2);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc]
//                == 0.3 + (0.6 + 0.1) * p_detach);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.6 * (1 - p_detach));  // loc is {0, 1}
//     loc[1] = 2;
//     BOOST_TEST(psv2.tensor[loc] == 0.1 * (1 - p_detach));  // loc is {0, 2}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.4;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    psv.tensor[loc] = 0.5;  // loc is {2, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    dt.forward(&edmans, &psv);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    loc[0] = 0;
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * (1 - p_detach));  // loc is {0, 1}
    loc[0] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.4 * (1 - p_detach));  // loc is {1, 1}
    loc[0] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.6 * (1 - p_detach));  // loc is {2, 1}
    // Distribution between empty states is of no importance. We only care about
    // the sum.
    double sum_empties = 0.0;
    loc[0] = 0;
    loc[1] = 0;
    sum_empties += psv.tensor[loc];  // loc is {0, 0}
    loc[0] = 1;
    sum_empties += psv.tensor[loc];  // loc is {1, 0};
    loc[0] = 2;
    sum_empties += psv.tensor[loc];  // loc is {2, 0};
    BOOST_TEST(sum_empties == 0.1 + 0.3 + 0.5 + (0.2 + 0.4 + 0.6) * p_detach);
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_multiple_edmans_test, *tolerance(TOL)) {
//     double p_detach = 0.05;
//     DetachTransition dt(p_detach);
//     int order = 2;
//     int* shape = new int[order];
//     shape[0] = 3;
//     shape[1] = 2;
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.1;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.2;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {1, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.4;  // loc is {1, 1}
//     loc[0] = 2;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.5;  // loc is {2, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.6;  // loc is {2, 1}
//     int edmans = 2;
//     dt.forward(tsr1, &edmans, &psv2);
//     // Just testing the ones with at least one lit amino acid here. See below
//     // for other tests.
//     loc[0] = 0;
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 * (1 - p_detach));  // loc is {0, 1}
//     loc[0] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.4 * (1 - p_detach));  // loc is {1, 1}
//     loc[0] = 2;
//     BOOST_TEST(psv2.tensor[loc] == 0.6 * (1 - p_detach));  // loc is {2, 1}
//     // Distribution between empty states is of no importance. We only care
//     about
//     // the sum.
//     double sum_empties = 0.0;
//     loc[0] = 0;
//     loc[1] = 0;
//     sum_empties += psv2.tensor[loc];  // loc is {0, 0}
//     loc[0] = 1;
//     sum_empties += psv2.tensor[loc];  // loc is {1, 0};
//     loc[0] = 2;
//     sum_empties += psv2.tensor[loc];  // loc is {2, 0};
//     BOOST_TEST(sum_empties == 0.1 + 0.3 + 0.5 + (0.2 + 0.4 + 0.6) *
//     p_detach); delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
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
    int edmans = 0;
    dt.forward(&edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {0, 0, 0}
    BOOST_TEST(psv.tensor[loc] == 0.1 + (0.2 + 0.3 + 0.4) * p_detach);
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 * (1 - p_detach));  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * (1 - p_detach));  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.4 * (1 - p_detach));  // loc is {0, 1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_multiple_dye_colors_test,
//                      *tolerance(TOL)) {
//     double p_detach = 0.05;
//     DetachTransition dt(p_detach);
//     int order = 3;
//     int* shape = new int[order];
//     shape[0] = 1;
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
//     int edmans = 0;
//     dt.forward(tsr1, &edmans, &psv2);
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     // loc is {0, 0, 0}
//     BOOST_TEST(psv2.tensor[loc] == 0.1 + (0.2 + 0.3 + 0.4) * p_detach);
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc]
//                == 0.2 * (1 - p_detach));  // loc is {0, 0, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     BOOST_TEST(psv2.tensor[loc]
//                == 0.3 * (1 - p_detach));  // loc is {0, 1, 0}
//     loc[2] = 1;
//     BOOST_TEST(psv2.tensor[loc]
//                == 0.4 * (1 - p_detach));  // loc is {0, 1, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(backward_in_place_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 1.0;  // loc is {0, 0}
    int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    BOOST_TEST(psv.tensor[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 1.0;  // loc is {0, 0}
    int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    BOOST_TEST(psv2.tensor[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {0, 1}
    int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.3 + (1 - p_detach) * 0.7);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
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
    int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.3 + (1 - p_detach) * 0.7);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.6;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.1;  // loc is {0, 2}
    int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.3 + (1 - p_detach) * 0.6);
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.3 + (1 - p_detach) * 0.1);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.6;  // loc is {0, 1}
    loc[1] = 2;
    psv1.tensor[loc] = 0.1;  // loc is {0, 2}
    int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.3 + (1 - p_detach) * 0.6);
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.3 + (1 - p_detach) * 0.1);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.88;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 0.88;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.4;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    psv.tensor[loc] = 0.88;  // loc is {2, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    dt.backward(psv, &edmans, &psv);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.88);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.88 + (1 - p_detach) * 0.2);
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.88);  // loc is {1, 0}
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.88 + (1 - p_detach) * 0.4);
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.88);  // loc is {2, 0}
    loc[1] = 1;
    // loc is {2, 1}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.88 + (1 - p_detach) * 0.6);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv1.tensor[loc] = 0.88;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.2;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv1.tensor[loc] = 0.88;  // loc is {1, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.4;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    psv1.tensor[loc] = 0.88;  // loc is {2, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    dt.backward(psv1, &edmans, &psv2);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.88);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.88 + (1 - p_detach) * 0.2);
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.88);  // loc is {1, 0}
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.88 + (1 - p_detach) * 0.4);
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.88);  // loc is {2, 0}
    loc[1] = 1;
    // loc is {2, 1}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.88 + (1 - p_detach) * 0.6);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
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
    int edmans = 0;
    dt.backward(psv, &edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1);  // loc is {0, 0, 0}
    loc[2] = 1;
    // loc is {0, 0, 1}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.1 + (1 - p_detach) * 0.2);
    loc[1] = 1;
    loc[2] = 0;
    // loc is {0, 1, 0}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.1 + (1 - p_detach) * 0.3);
    loc[2] = 1;
    // loc is {0, 1, 1}
    BOOST_TEST(psv.tensor[loc] == p_detach * 0.1 + (1 - p_detach) * 0.4);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_multiple_dye_colors_test,
                     *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
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
    int edmans = 0;
    dt.backward(psv1, &edmans, &psv2);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.1);  // loc is {0, 0, 0}
    loc[2] = 1;
    // loc is {0, 0, 1}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.1 + (1 - p_detach) * 0.2);
    loc[1] = 1;
    loc[2] = 0;
    // loc is {0, 1, 0}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.1 + (1 - p_detach) * 0.3);
    loc[2] = 1;
    // loc is {0, 1, 1}
    BOOST_TEST(psv2.tensor[loc] == p_detach * 0.1 + (1 - p_detach) * 0.4);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(improve_fit_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector fpsv(order, shape);
    PeptideStateVector bpsv(order, shape);
    PeptideStateVector nbpsv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    fpsv.tensor[loc] = 0.31;  // loc is {0, 0}
    bpsv.tensor[loc] = 0.32;
    nbpsv.tensor[loc] = 0.33;
    loc[1] = 1;
    fpsv.tensor[loc] = 0.71;  // loc is {0, 1}
    bpsv.tensor[loc] = 0.72;
    nbpsv.tensor[loc] = 0.73;
    loc[1] = 2;
    fpsv.tensor[loc] = 0.91;  // loc is {0, 2}
    bpsv.tensor[loc] = 0.92;
    nbpsv.tensor[loc] = 0.93;
    delete[] loc;
    int edmans = 0;
    double probability = 0.31 * 0.32 + 0.71 * 0.72 + 0.91 * 0.92;
    ErrorModelFitter emf;
    dt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &emf);
    BOOST_TEST(emf.p_detach_fit.get()
               == (0.71 * p_detach * 0.33 + 0.91 * p_detach * 0.33)
                          / (0.71 * 0.72 + 0.91 * 0.92));
}

BOOST_AUTO_TEST_CASE(improve_fit_twice_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
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
    bpsv1.tensor[loc] = 0.32;
    nbpsv1.tensor[loc] = 0.33;
    fpsv2.tensor[loc] = 0.231;
    bpsv2.tensor[loc] = 0.232;
    nbpsv2.tensor[loc] = 0.233;
    loc[1] = 1;
    fpsv1.tensor[loc] = 0.71;  // loc is {0, 1}
    bpsv1.tensor[loc] = 0.72;
    nbpsv1.tensor[loc] = 0.73;
    fpsv2.tensor[loc] = 0.271;
    bpsv2.tensor[loc] = 0.272;
    nbpsv2.tensor[loc] = 0.273;
    loc[1] = 2;
    fpsv1.tensor[loc] = 0.91;  // loc is {0, 2}
    bpsv1.tensor[loc] = 0.92;
    nbpsv1.tensor[loc] = 0.93;
    fpsv2.tensor[loc] = 0.291;
    bpsv2.tensor[loc] = 0.292;
    nbpsv2.tensor[loc] = 0.293;
    delete[] loc;
    int edmans = 0;
    double prob1 = 0.31 * 0.32 + 0.71 * 0.72 + 0.91 * 0.92;
    double prob2 = 0.231 * 0.232 + 0.271 * 0.272 + 0.291 * 0.292;
    ErrorModelFitter emf;
    dt.improve_fit(fpsv1, bpsv1, nbpsv1, edmans, prob1, &emf);
    dt.improve_fit(fpsv2, bpsv2, nbpsv2, edmans, prob2, &emf);
    BOOST_TEST(
            emf.p_detach_fit.get()
            == ((0.71 * p_detach * 0.33 + 0.91 * p_detach * 0.33) / prob1
                + (0.271 * p_detach * 0.233 + 0.291 * p_detach * 0.233) / prob2)
                       / ((0.71 * 0.72 + 0.91 * 0.92) / prob1
                          + (0.271 * 0.272 + 0.291 * 0.292) / prob2));
}

BOOST_AUTO_TEST_SUITE_END()  // detach_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
