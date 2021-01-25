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
#include "hmm.h"

// Standard C++ library headers:
#include <cmath>
#include <functional>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/binomial-transition.h"
#include "hmm/detach-transition.h"
#include "hmm/dye-seq-precomputations.h"
#include "hmm/edman-transition.h"
#include "hmm/emission.h"
#include "hmm/radiometry-precomputations.h"
#include "hmm/start.h"
#include "hmm/universal-precomputations.h"
#include "tensor/tensor.h"

namespace fluoroseq {

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

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double p_edman_failure = 0.01;
    double p_detach = 0.02;
    double p_bleach = 0.03;
    double p_dud = 0.04;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = 1.0;
    double sigma = 0.05;
    ErrorModel em(
            p_edman_failure, p_detach, p_bleach, p_dud, dist_type, mu, sigma);
    int num_channels = 2;
    int max_num_dyes = 3;
    UniversalPrecomputations universal_precomputations(em, num_channels);
    universal_precomputations.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 4;
    DyeSeq ds(num_channels, ".1.0.1.0.1");  // two in ch 0, three in ch 1.
    DyeSeqPrecomputations dye_seq_precomputations(
            ds, em, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 1.0;
    r(0, 1) = 1.0;
    r(1, 0) = 1.0;
    r(1, 1) = 1.0;
    r(2, 0) = 1.0;
    r(2, 1) = 1.0;
    r(3, 0) = 1.0;
    r(3, 1) = 1.0;
    RadiometryPrecomputations radiometry_precomputations(r, em, max_num_dyes);
    HMM hmm(num_timesteps,
            num_channels,
            dye_seq_precomputations,
            radiometry_precomputations,
            universal_precomputations);
    BOOST_ASSERT(hmm.tensor_shape.size() == 1 + num_channels);
    BOOST_TEST(hmm.tensor_shape[0] == num_timesteps);
    BOOST_TEST(hmm.tensor_shape[1] == 2 + 1);  // extra is to have 0 & num dyes.
    BOOST_TEST(hmm.tensor_shape[2] == 3 + 1);  // extra is to have 0 & num dyes.
    BOOST_ASSERT(hmm.steps.size()
                 == (3 + num_channels) * (num_timesteps - 1) + 2
                            + num_channels);
    vector<const Step*>::iterator step = hmm.steps.begin();
    BOOST_TEST(*step == &universal_precomputations.start);
    step++;
    BOOST_TEST(*step == &universal_precomputations.dud_transitions[0]);
    step++;
    BOOST_TEST(*step == &universal_precomputations.dud_transitions[1]);
    step++;
    BOOST_TEST(*step == &radiometry_precomputations.emission);
    step++;
    BOOST_TEST(*step == &universal_precomputations.detach_transition);
    step++;
    BOOST_TEST(*step == &universal_precomputations.bleach_transitions[0]);
    step++;
    BOOST_TEST(*step == &universal_precomputations.bleach_transitions[1]);
    step++;
    BOOST_TEST(*step == &dye_seq_precomputations.edman_transition);
    step++;
    BOOST_TEST(*step == &radiometry_precomputations.emission);
    step++;
    BOOST_TEST(*step == &universal_precomputations.detach_transition);
    step++;
    BOOST_TEST(*step == &universal_precomputations.bleach_transitions[0]);
    step++;
    BOOST_TEST(*step == &universal_precomputations.bleach_transitions[1]);
    step++;
    BOOST_TEST(*step == &dye_seq_precomputations.edman_transition);
    step++;
    BOOST_TEST(*step == &radiometry_precomputations.emission);
}

BOOST_AUTO_TEST_CASE(probability_trivial_test, *tolerance(TOL)) {
    double p_edman_failure = 0.0;
    double p_detach = 0.0;
    double p_bleach = 0.0;
    double p_dud = 0.0;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = 0.05;
    ErrorModel em(
            p_edman_failure, p_detach, p_bleach, p_dud, dist_type, mu, sigma);
    int num_channels = 0;
    int max_num_dyes = 0;
    UniversalPrecomputations up(em, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 1;
    DyeSeq ds(num_channels, "");  // 0 in all channels.
    DyeSeqPrecomputations dsp(ds, em, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    RadiometryPrecomputations rp(r, em, max_num_dyes);
    HMM hmm(num_timesteps, num_channels, dsp, rp, up);
    BOOST_TEST(hmm.probability() == 1.0);
}

BOOST_AUTO_TEST_CASE(probability_sum_to_one_test, *tolerance(TOL)) {
    // The idea here is that if the emission just causes a multiplication by one
    // no matter what, then the transitions shouldn't matter; they should just
    // be reallocating a probability of one between different states.
    double p_edman_failure = 0.45;
    double p_detach = 0.15;
    double p_bleach = 0.35;
    double p_dud = 0.25;
    DistributionType dist_type = DistributionType::OVERRIDE;
    double mu = 1.0;
    double sigma = 0.05;
    ErrorModel em(
            p_edman_failure, p_detach, p_bleach, p_dud, dist_type, mu, sigma);
    int num_channels = 2;
    int max_num_dyes = 5;
    UniversalPrecomputations up(em, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 3;
    DyeSeq ds(num_channels, "10.01111");  // two in ch 0, five in ch 1.
    DyeSeqPrecomputations dsp(ds, em, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 0.0;
    r(0, 1) = 0.1;
    r(1, 0) = 1.0;
    r(1, 1) = 1.1;
    r(2, 0) = 2.0;
    r(2, 1) = 2.1;
    RadiometryPrecomputations rp(r, em, max_num_dyes);
    HMM hmm(num_timesteps, num_channels, dsp, rp, up);
    BOOST_TEST(hmm.probability() == 1.0);
}

BOOST_AUTO_TEST_CASE(probability_more_involved_test, *tolerance(TOL)) {
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on November 10, 2020.
    double p_edman_failure = 0.06;
    double p_detach = 0.05;
    double p_bleach = 0.05;
    double p_dud = 0.07;
    DistributionType dist_type = DistributionType::LOGNORMAL;
    double mu = 1.0;
    double sigma = 0.16;
    ErrorModel em(
            p_edman_failure, p_detach, p_bleach, p_dud, dist_type, mu, sigma);
    int num_channels = 2;
    int max_num_dyes = 5;
    UniversalPrecomputations up(em, num_channels);
    up.set_max_num_dyes(max_num_dyes);
    int num_timesteps = 3;
    DyeSeq ds(num_channels, "10.01111");  // two in ch 0, five in ch 1.
    DyeSeqPrecomputations dsp(ds, em, num_timesteps, num_channels);
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 5.0;
    r(0, 1) = 2.0;
    r(1, 0) = 5.0;
    r(1, 1) = 1.0;
    r(2, 0) = 4.0;
    r(2, 1) = 1.0;
    RadiometryPrecomputations rp(r, em, max_num_dyes);
    HMM hmm(num_timesteps, num_channels, dsp, rp, up);
    BOOST_TEST(hmm.probability() == 3.2324422559808915e-23);
}

BOOST_AUTO_TEST_SUITE_END()  // hmm_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace fluoroseq
