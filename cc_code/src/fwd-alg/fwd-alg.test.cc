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
#include "fwd-alg.h"

// Standard C++ library headers:
#include <cmath>
#include <functional>

// Local project headers:
#include "fwd-alg/binomial-transition.h"
#include "fwd-alg/detach-transition.h"
#include "fwd-alg/edman-transition.h"
#include "fwd-alg/emission.h"
#include "fwd-alg/initialization.h"
#include "fwd-alg/summation.h"
#include "tensor/tensor.h"

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
using std::exp;
using std::function;
using std::log;
using std::sqrt;
const double PI = 3.141592653589793238;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(fwd_alg_suite)
BOOST_AUTO_TEST_SUITE(fwd_alg_suite)

BOOST_AUTO_TEST_CASE(trivial_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 0;
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    Tensor tsr(order, shape);
    delete[] shape;
    Initialization init;
    Radiometry rad(num_timesteps, num_channels);
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0;
    };
    Emission emit(rad, max_num_dyes, pdf);
    double p_detach = 0.0;
    DetachTransition detach(p_detach);
    double p_dud = 0.0;
    BinomialTransition dud(p_dud);
    double p_bleach = 0.0;
    BinomialTransition bleach(p_bleach);
    double p_edman_failure = 0.0;
    DyeSeq dye_seq(num_channels, "");
    DyeTrack dye_track(num_timesteps, num_channels, dye_seq);
    EdmanTransition edman(p_edman_failure, dye_seq, dye_track);
    Summation sum;
    BOOST_TEST(fwd_alg(&tsr,
                       num_timesteps,
                       num_channels,
                       init,
                       emit,
                       detach,
                       dud,
                       bleach,
                       edman,
                       sum)
               == 1.0);
}

BOOST_AUTO_TEST_CASE(sum_to_one_test, *tolerance(TOL)) {
    // The idea here is that if the emission just causes a multiplication by one
    // no matter what, then the transitions shouldn't matter; they should just
    // be reallocating a probability of one between different states.
    int num_timesteps = 3;
    int num_channels = 2;
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = 3;  // num 0s in dye_seq below plus one extra.
    shape[2] = 6;  // num 1s in dye_seq below plus one extra.
    Tensor tsr(order, shape);
    delete[] shape;
    Initialization init;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    rad(2, 0) = 2.0;
    rad(2, 1) = 2.1;
    int max_num_dyes = 5;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0;
    };
    Emission emit(rad, max_num_dyes, pdf);
    double p_detach = 0.15;
    DetachTransition detach(p_detach);
    double p_dud = 0.25;
    BinomialTransition dud(p_dud);
    dud.reserve(5);
    double p_bleach = 0.35;
    BinomialTransition bleach(p_bleach);
    bleach.reserve(5);
    double p_edman_failure = 0.45;
    DyeSeq dye_seq(num_channels, "10.01111");
    DyeTrack dye_track(num_timesteps, num_channels, dye_seq);
    EdmanTransition edman(p_edman_failure, dye_seq, dye_track);
    Summation sum;
    BOOST_TEST(fwd_alg(&tsr,
                       num_timesteps,
                       num_channels,
                       init,
                       emit,
                       detach,
                       dud,
                       bleach,
                       edman,
                       sum)
               == 1.0);
}

BOOST_AUTO_TEST_CASE(more_involved_test, *tolerance(TOL)) {
    // This is essentially a "no change" test. It assumes that the function was
    // giving the correct result on November 10, 2020.
    int num_timesteps = 3;
    int num_channels = 2;
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = 3;  // num 0s in dye_seq below plus one extra.
    shape[2] = 6;  // num 1s in dye_seq below plus one extra.
    Tensor tsr(order, shape);
    delete[] shape;
    Initialization init;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 5.0;
    rad(0, 1) = 2.0;
    rad(1, 0) = 5.0;
    rad(1, 1) = 1.0;
    rad(2, 0) = 4.0;
    rad(2, 1) = 1.0;
    int max_num_dyes = 5;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        double sigma = 0.16;
        if (state > 0) {
            if (observed == 0.0) {
                return 0.0;
            } else {
                double unit_obs = observed;
                double offset = log(unit_obs) - log((double)state);
                return (1.0 / (sigma * sqrt(2.0 * PI)) / unit_obs)
                       * exp(-(offset * offset) / (2.0 * sigma * sigma));
            }
        } else {
            if (observed == 0.0) {
                return 1.0;
            } else {
                return 0.0;
            }
        }
    };
    Emission emit(rad, max_num_dyes, pdf);
    double p_detach = 0.05;
    DetachTransition detach(p_detach);
    double p_dud = 0.07;
    BinomialTransition dud(p_dud);
    dud.reserve(5);
    double p_bleach = 0.05;
    BinomialTransition bleach(p_bleach);
    bleach.reserve(5);
    double p_edman_failure = 0.06;
    DyeSeq dye_seq(num_channels, "10.01111");
    DyeTrack dye_track(num_timesteps, num_channels, dye_seq);
    EdmanTransition edman(p_edman_failure, dye_seq, dye_track);
    Summation sum;
    BOOST_TEST(fwd_alg(&tsr,
                       num_timesteps,
                       num_channels,
                       init,
                       emit,
                       detach,
                       dud,
                       bleach,
                       edman,
                       sum)
               == 3.2324422559808915e-23);
}

BOOST_AUTO_TEST_SUITE_END()  // fwd_alg_suite
BOOST_AUTO_TEST_SUITE_END()  // fwd_alg_suite

}  // namespace fluoroseq
