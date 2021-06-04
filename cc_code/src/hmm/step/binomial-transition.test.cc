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
#include "binomial-transition.h"

// Local project headers:
#include "parameterization/fit/parameter-fitter.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

// BinomialTransition is abstract, so we need to override undefined methods.
class TestableBinomialTransition : public BinomialTransition {
public:
    TestableBinomialTransition(double q, int channel)
            : BinomialTransition(q, channel) {}
    using BinomialTransition::improve_fit;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
};

void TestableBinomialTransition::improve_fit(
        const PeptideStateVector& forward_psv,
        const PeptideStateVector& backward_psv,
        const PeptideStateVector& next_backward_psv,
        int num_edmans,
        double probability,
        SequencingModelFitter* fitter) const {}

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(binomial_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double q = 0.2;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    BOOST_TEST(bt.q == q);
    BOOST_TEST(bt.channel == channel);
}

BOOST_AUTO_TEST_CASE(reserve_zero_test, *tolerance(TOL)) {
    double q = 0.2;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
}

BOOST_AUTO_TEST_CASE(reserve_one_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
}

BOOST_AUTO_TEST_CASE(reserve_two_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
    BOOST_TEST(bt.prob(2, 0) == q * q);
    BOOST_TEST(bt.prob(2, 1) == 2 * q * p);
    BOOST_TEST(bt.prob(2, 2) == p * p);
}

BOOST_AUTO_TEST_CASE(reserve_three_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(3);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
    BOOST_TEST(bt.prob(2, 0) == q * q);
    BOOST_TEST(bt.prob(2, 1) == 2 * q * p);
    BOOST_TEST(bt.prob(2, 2) == p * p);
    BOOST_TEST(bt.prob(3, 0) == q * q * q);
    BOOST_TEST(bt.prob(3, 1) == 3 * q * q * p);
    BOOST_TEST(bt.prob(3, 2) == 3 * q * p * p);
    BOOST_TEST(bt.prob(3, 3) == p * p * p);
}

BOOST_AUTO_TEST_CASE(reserve_no_shrink_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(3);
    bt.reserve(1);  // This operation should be silently ignored.
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
    BOOST_TEST(bt.prob(2, 0) == q * q);
    BOOST_TEST(bt.prob(2, 1) == 2 * q * p);
    BOOST_TEST(bt.prob(2, 2) == p * p);
    BOOST_TEST(bt.prob(3, 0) == q * q * q);
    BOOST_TEST(bt.prob(3, 1) == 3 * q * q * p);
    BOOST_TEST(bt.prob(3, 2) == 3 * q * p * p);
    BOOST_TEST(bt.prob(3, 3) == p * p * p);
}

BOOST_AUTO_TEST_CASE(forward_in_place_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
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
    bt.forward(&edmans, &psv);
    BOOST_TEST(psv.tensor[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_trivial_test, *tolerance(TOL)) {
//     double q = 0.05;
//     double p = 0.95;
//     int channel = 0;
//     TestableBinomialTransition bt(q, channel);
//     bt.reserve(0);
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
//     bt.forward(tsr1, &edmans, &psv2);
//     BOOST_TEST(psv2.tensor[loc] == 1.0);  // loc is {0, 0}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
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
    bt.forward(&edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 + 0.7 * q);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p);  // loc is {0, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_basic_transition_test, *tolerance(TOL))
// {
//     double q = 0.05;
//     double p = 0.95;
//     int channel = 0;
//     TestableBinomialTransition bt(q, channel);
//     bt.reserve(1);
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
//     bt.forward(tsr1, &edmans, &psv2);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 + 0.7 * q);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p);  // loc is {0, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_bigger_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.3;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.7;  // loc is {0, 2}
    int edmans = 0;
    bt.forward(&edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    // loc is {0, 0}
    BOOST_TEST(psv.tensor[loc] == 0.2 + 0.3 * q + 0.7 * q * q);
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p + 0.7 * 2 * q * p);  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p * p);  // loc is {0, 2}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_bigger_transition_test, *tolerance(TOL))
// {
//     double q = 0.05;
//     double p = 0.95;
//     int channel = 0;
//     TestableBinomialTransition bt(q, channel);
//     bt.reserve(2);
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
//     psv1.tensor[loc] = 0.2;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 1}
//     loc[1] = 2;
//     psv1.tensor[loc] = 0.7;  // loc is {0, 2}
//     int edmans = 0;
//     bt.forward(tsr1, &edmans, &psv2);
//     loc[0] = 0;
//     loc[1] = 0;
//     // loc is {0, 0}
//     BOOST_TEST(psv2.tensor[loc] == 0.2 + 0.3 * q + 0.7 * q * q);
//     loc[1] = 1;
//     // loc is {0, 1}
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p + 0.7 * 2 * q * p);
//     loc[1] = 2;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p * p);  // loc is {0, 2}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.8;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    psv.tensor[loc] = 0.4;  // loc is {2, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    bt.forward(&edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.2 + 0.8 * q);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.8 * p);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 + 0.7 * q);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p);  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.4 + 0.6 * q);  // loc is {2, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.6 * p);  // loc is {2, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_multiple_edmans_test, *tolerance(TOL)) {
//     double q = 0.05;
//     double p = 0.95;
//     int channel = 0;
//     TestableBinomialTransition bt(q, channel);
//     bt.reserve(2);
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
//     psv1.tensor[loc] = 0.2;  // loc is {0, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.8;  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {1, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.7;  // loc is {1, 1}
//     loc[0] = 2;
//     loc[1] = 0;
//     psv1.tensor[loc] = 0.4;  // loc is {2, 0}
//     loc[1] = 1;
//     psv1.tensor[loc] = 0.6;  // loc is {2, 1}
//     int edmans = 2;
//     bt.forward(tsr1, &edmans, &psv2);
//     loc[0] = 0;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 + 0.8 * q);  // loc is {0, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.8 * p);  // loc is {0, 1}
//     loc[0] = 1;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 + 0.7 * q);  // loc is {1, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p);  // loc is {1, 1}
//     loc[0] = 2;
//     loc[1] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.4 + 0.6 * q);  // loc is {2, 0}
//     loc[1] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.6 * p);  // loc is {2, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(forward_in_place_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 4;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.4;  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    psv.tensor[loc] = 0.5;  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.6;  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    psv.tensor[loc] = 0.7;  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.8;  // loc is {0, 1, 1, 1}
    int edmans = 0;
    bt.forward(&edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1 + 0.3 * q);  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2 + 0.4 * q);  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3 * p);  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.4 * p);  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.5 + 0.7 * q);  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.6 + 0.8 * q);  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.7 * p);  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.8 * p);  // loc is {0, 1, 1, 1}
    delete[] loc;
}

// BOOST_AUTO_TEST_CASE(forward_new_tsr_other_dye_colors_test, *tolerance(TOL))
// {
//     double q = 0.05;
//     double p = 0.95;
//     int channel = 1;  // corresponds to 2nd dim of tensor
//     TestableBinomialTransition bt(q, channel);
//     bt.reserve(1);
//     int order = 4;
//     int* shape = new int[order];
//     shape[0] = 1;
//     shape[1] = 2;
//     shape[2] = 2;
//     shape[3] = 2;
//     PeptideStateVector psv1(order, shape);
//     PeptideStateVector psv2(order, shape);
//     delete[] shape;
//     int* loc = new int[order];
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     loc[3] = 0;
//     psv1.tensor[loc] = 0.1;  // loc is {0, 0, 0, 0}
//     loc[3] = 1;
//     psv1.tensor[loc] = 0.2;  // loc is {0, 0, 0, 1}
//     loc[2] = 1;
//     loc[3] = 0;
//     psv1.tensor[loc] = 0.3;  // loc is {0, 0, 1, 0}
//     loc[3] = 1;
//     psv1.tensor[loc] = 0.4;  // loc is {0, 0, 1, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     loc[3] = 0;
//     psv1.tensor[loc] = 0.5;  // loc is {0, 1, 0, 0}
//     loc[3] = 1;
//     psv1.tensor[loc] = 0.6;  // loc is {0, 1, 0, 1}
//     loc[2] = 1;
//     loc[3] = 0;
//     psv1.tensor[loc] = 0.7;  // loc is {0, 1, 1, 0}
//     loc[3] = 1;
//     psv1.tensor[loc] = 0.8;  // loc is {0, 1, 1, 1}
//     int edmans = 0;
//     bt.forward(tsr1, &edmans, &psv2);
//     loc[0] = 0;
//     loc[1] = 0;
//     loc[2] = 0;
//     loc[3] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.1 + 0.3 * q);  // loc is {0, 0, 0, 0}
//     loc[3] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.2 + 0.4 * q);  // loc is {0, 0, 0, 1}
//     loc[2] = 1;
//     loc[3] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.3 * p);  // loc is {0, 0, 1, 0}
//     loc[3] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.4 * p);  // loc is {0, 0, 1, 1}
//     loc[1] = 1;
//     loc[2] = 0;
//     loc[3] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.5 + 0.7 * q);  // loc is {0, 1, 0, 0}
//     loc[3] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.6 + 0.8 * q);  // loc is {0, 1, 0, 1}
//     loc[2] = 1;
//     loc[3] = 0;
//     BOOST_TEST(psv2.tensor[loc] == 0.7 * p);  // loc is {0, 1, 1, 0}
//     loc[3] = 1;
//     BOOST_TEST(psv2.tensor[loc] == 0.8 * p);  // loc is {0, 1, 1, 1}
//     delete[] loc;
// }

BOOST_AUTO_TEST_CASE(backward_in_place_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
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
    bt.backward(psv, &edmans, &psv);
    BOOST_TEST(psv.tensor[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
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
    bt.backward(psv1, &edmans, &psv2);
    BOOST_TEST(psv2.tensor[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
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
    bt.backward(psv, &edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == q * 0.3 + p * 0.7);  // loc is {0, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
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
    bt.backward(psv1, &edmans, &psv2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc] == q * 0.3 + p * 0.7);  // loc is {0, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_bigger_transition_test,
                     *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.3;  // loc is {0, 1}
    loc[1] = 2;
    psv.tensor[loc] = 0.7;  // loc is {0, 2}
    int edmans = 0;
    bt.backward(psv, &edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.2);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == q * 0.2 + p * 0.3);  // loc is {0, 1}
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(psv.tensor[loc] == q * q * 0.2 + 2 * q * p * 0.3 + p * p * 0.7);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_bigger_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
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
    psv1.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.3;  // loc is {0, 1}
    loc[1] = 2;
    psv1.tensor[loc] = 0.7;  // loc is {0, 2}
    int edmans = 0;
    bt.backward(psv1, &edmans, &psv2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.2);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc] == q * 0.2 + p * 0.3);  // loc is {0, 1}
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(psv2.tensor[loc] == q * q * 0.2 + 2 * q * p * 0.3 + p * p * 0.7);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    psv.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.8;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv.tensor[loc] = 0.3;  // loc is {1, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.7;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    psv.tensor[loc] = 0.4;  // loc is {2, 0}
    loc[1] = 1;
    psv.tensor[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    bt.backward(psv, &edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.2);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == q * 0.2 + p * 0.8);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.3);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == q * 0.3 + p * 0.7);  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.4);  // loc is {2, 0}
    loc[1] = 1;
    BOOST_TEST(psv.tensor[loc] == q * 0.4 + p * 0.6);  // loc is {2, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
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
    psv1.tensor[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.8;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {1, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.7;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    psv1.tensor[loc] = 0.4;  // loc is {2, 0}
    loc[1] = 1;
    psv1.tensor[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    bt.backward(psv1, &edmans, &psv2);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.2);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc] == q * 0.2 + p * 0.8);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.3);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc] == q * 0.3 + p * 0.7);  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.4);  // loc is {2, 0}
    loc[1] = 1;
    BOOST_TEST(psv2.tensor[loc] == q * 0.4 + p * 0.6);  // loc is {2, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_in_place_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 4;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    psv.tensor[loc] = 0.1;  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.2;  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    psv.tensor[loc] = 0.3;  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.4;  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    psv.tensor[loc] = 0.5;  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.6;  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    psv.tensor[loc] = 0.7;  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    psv.tensor[loc] = 0.8;  // loc is {0, 1, 1, 1}
    int edmans = 0;
    bt.backward(psv, &edmans, &psv);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.1);  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.2);  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == q * 0.1 + p * 0.3);  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == q * 0.2 + p * 0.4);  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == 0.5);  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == 0.6);  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(psv.tensor[loc] == q * 0.5 + p * 0.7);  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    BOOST_TEST(psv.tensor[loc] == q * 0.6 + p * 0.8);  // loc is {0, 1, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 4;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    psv1.tensor[loc] = 0.1;  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    psv1.tensor[loc] = 0.2;  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    psv1.tensor[loc] = 0.3;  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    psv1.tensor[loc] = 0.4;  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    psv1.tensor[loc] = 0.5;  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    psv1.tensor[loc] = 0.6;  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    psv1.tensor[loc] = 0.7;  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    psv1.tensor[loc] = 0.8;  // loc is {0, 1, 1, 1}
    int edmans = 0;
    bt.backward(psv1, &edmans, &psv2);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.1);  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    BOOST_TEST(psv2.tensor[loc] == 0.2);  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(psv2.tensor[loc] == q * 0.1 + p * 0.3);  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    BOOST_TEST(psv2.tensor[loc] == q * 0.2 + p * 0.4);  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(psv2.tensor[loc] == 0.5);  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    BOOST_TEST(psv2.tensor[loc] == 0.6);  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(psv2.tensor[loc] == q * 0.5 + p * 0.7);  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    BOOST_TEST(psv2.tensor[loc] == q * 0.6 + p * 0.8);  // loc is {0, 1, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(improve_fit_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector fpsv(order, shape);
    PeptideStateVector bpsv(order, shape);
    PeptideStateVector nbpsv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    fpsv.tensor[loc] = 1.0;  // loc is {0, 0}
    bpsv.tensor[loc] = 1.0;
    nbpsv.tensor[loc] = 1.0;
    delete[] loc;
    int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.numerator == 0.0);
    BOOST_TEST(pf.denominator == 0.0);
}

BOOST_AUTO_TEST_CASE(improve_fit_basic_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
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
    delete[] loc;
    int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get() == (0.71 * q * 0.33) / (0.71 * 0.72));
}

BOOST_AUTO_TEST_CASE(improve_fit_bigger_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
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
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 + 0.91 * (q * p * 2.0) * 0.73
                   + 0.91 * (q * q) * 0.33 * 2.0)
                          / (0.71 * 0.72 + 0.91 * 0.92 * 2.0));
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
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
    fpsv.tensor[loc] = 0.31;  // loc is {0, 0}
    bpsv.tensor[loc] = 0.32;
    nbpsv.tensor[loc] = 0.33;
    loc[1] = 1;
    fpsv.tensor[loc] = 0.71;  // loc is {0, 1}
    bpsv.tensor[loc] = 0.72;
    nbpsv.tensor[loc] = 0.73;
    loc[0] = 1;
    loc[1] = 0;
    fpsv.tensor[loc] = 0.41;  // loc is {1, 0}
    bpsv.tensor[loc] = 0.42;
    nbpsv.tensor[loc] = 0.43;
    loc[1] = 1;
    fpsv.tensor[loc] = 0.81;  // loc is {1, 1}
    bpsv.tensor[loc] = 0.82;
    nbpsv.tensor[loc] = 0.83;
    delete[] loc;
    int edmans = 1;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 + 0.81 * q * 0.43)
                          / (0.71 * 0.72 + 0.81 * 0.82));
}

BOOST_AUTO_TEST_CASE(improve_fit_other_dye_color_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector fpsv(order, shape);
    PeptideStateVector bpsv(order, shape);
    PeptideStateVector nbpsv(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    fpsv.tensor[loc] = 0.31;  // loc is {0, 0, 0}
    bpsv.tensor[loc] = 0.32;
    nbpsv.tensor[loc] = 0.33;
    loc[2] = 1;
    fpsv.tensor[loc] = 0.231;  // loc is {0, 0, 1}
    bpsv.tensor[loc] = 0.232;
    nbpsv.tensor[loc] = 0.233;
    loc[1] = 1;
    loc[2] = 0;
    fpsv.tensor[loc] = 0.71;  // loc is {0, 1, 0}
    bpsv.tensor[loc] = 0.72;
    nbpsv.tensor[loc] = 0.73;
    loc[2] = 1;
    fpsv.tensor[loc] = 0.271;  // loc is {0, 1, 1}
    bpsv.tensor[loc] = 0.272;
    nbpsv.tensor[loc] = 0.273;
    delete[] loc;
    int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 + 0.271 * q * 0.233)
                          / (0.71 * 0.72 + 0.271 * 0.272));
}

BOOST_AUTO_TEST_CASE(improve_fit_different_probability_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
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
    delete[] loc;
    int edmans = 0;
    double probability = 0.123456789;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get() == (0.71 * q * 0.33) / (0.71 * 0.72));
}

BOOST_AUTO_TEST_CASE(improve_fit_twice_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
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
    delete[] loc;
    int edmans = 0;
    double prob1 = 0.123456789;
    double prob2 = 0.987654321;
    ParameterFitter pf;
    bt.improve_fit(fpsv1, bpsv1, nbpsv1, edmans, prob1, &pf);
    bt.improve_fit(fpsv2, bpsv2, nbpsv2, edmans, prob2, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 / prob1 + 0.271 * q * 0.233 / prob2)
                          / (0.71 * 0.72 / prob1 + 0.271 * 0.272 / prob2));
}

BOOST_AUTO_TEST_SUITE_END()  // binomial_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
