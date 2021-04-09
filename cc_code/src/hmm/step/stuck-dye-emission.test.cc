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
#include "stuck-dye-emission.h"

// Standard C++ library headers:
#include <functional>

// External headers:
#include "fakeit.hpp"

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "test-util/fakeit.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using fakeit::Fake;
using fakeit::Mock;
using fakeit::Verify;
using fakeit::VerifyNoOtherInvocations;
using std::function;
using whatprot::test_util::Close;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(stuck_dye_emission_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 0.5;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.num_channels == num_channels);
    BOOST_TEST(e.channel == channel);
    BOOST_TEST(e.values[0] == 0.5);
}

BOOST_AUTO_TEST_CASE(prob_multiple_timesteps_test, *tolerance(TOL)) {
    int num_timesteps = 3;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 1));
    BOOST_TEST(e.prob(1, 0, 0) == pdf(1.0, 0));
    BOOST_TEST(e.prob(1, 0, 0) == pdf(1.0, 1));
    BOOST_TEST(e.prob(2, 0, 0) == pdf(2.0, 0));
    BOOST_TEST(e.prob(2, 0, 0) == pdf(2.0, 1));
}

BOOST_AUTO_TEST_CASE(prob_multiple_timesteps_const_test, *tolerance(TOL)) {
    int num_timesteps = 3;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    const StuckDyeEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 1));
    BOOST_TEST(ce.prob(1, 0, 0) == pdf(1.0, 0));
    BOOST_TEST(ce.prob(1, 0, 0) == pdf(1.0, 1));
    BOOST_TEST(ce.prob(2, 0, 0) == pdf(2.0, 0));
    BOOST_TEST(ce.prob(2, 0, 0) == pdf(2.0, 1));
}

BOOST_AUTO_TEST_CASE(prob_multiple_channels_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 1));
    BOOST_TEST(e.prob(0, 1, 0) == pdf(0.1, 0));
    BOOST_TEST(e.prob(0, 1, 0) == pdf(0.1, 1));
}

BOOST_AUTO_TEST_CASE(prob_multiple_channels_const_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    const StuckDyeEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 1));
    BOOST_TEST(ce.prob(0, 1, 0) == pdf(0.1, 0));
    BOOST_TEST(ce.prob(0, 1, 0) == pdf(0.1, 1));
}

BOOST_AUTO_TEST_CASE(forward_basic_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye == 0.3 * pdf(1.0, 0));
    BOOST_TEST(sdsv.dye == 0.7 * pdf(1.0, 1));
}

BOOST_AUTO_TEST_CASE(forward_multiple_timestamps_test, *tolerance(TOL)) {
    int num_timesteps = 2;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(1, 0) = 1.1;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 1;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye == 0.3 * pdf(1.1, 0));
    BOOST_TEST(sdsv.dye == 0.7 * pdf(1.1, 1));
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye == 0.3 * pdf(1.0, 0) * pdf(1.1, 0));
    BOOST_TEST(sdsv.dye == 0.7 * pdf(1.0, 1) * pdf(1.1, 0));
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_other_channel_test,
                     *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    int channel = 1;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye == 0.3 * pdf(1.0, 0) * pdf(1.1, 0));
    BOOST_TEST(sdsv.dye == 0.7 * pdf(1.0, 0) * pdf(1.1, 1));
}

BOOST_AUTO_TEST_CASE(backward_basic_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye == 0.3 * pdf(1.0, 0));
    BOOST_TEST(output.dye == 0.7 * pdf(1.0, 1));
}

BOOST_AUTO_TEST_CASE(backward_multiple_timestamps_test, *tolerance(TOL)) {
    int num_timesteps = 2;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(1, 0) = 1.1;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 1;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye == 0.3 * pdf(1.1, 0));
    BOOST_TEST(output.dye == 0.7 * pdf(1.1, 1));
}

BOOST_AUTO_TEST_CASE(backward_multiple_channels_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye == 0.3 * pdf(1.0, 0) * pdf(1.1, 0));
    BOOST_TEST(output.dye == 0.7 * pdf(1.0, 1) * pdf(1.1, 0));
}

BOOST_AUTO_TEST_CASE(backward_multiple_channels_other_channel_test,
                     *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    int channel = 1;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye == 0.3 * pdf(1.0, 0) * pdf(1.1, 0));
    BOOST_TEST(output.dye == 0.7 * pdf(1.0, 0) * pdf(1.1, 1));
}

BOOST_AUTO_TEST_CASE(improve_fit_basic_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.2;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    double probability = 0.98765;
    StuckDyeStateVector forward_sdsv;
    forward_sdsv.dye = 0.7;
    forward_sdsv.no_dye = 0.3;
    StuckDyeStateVector backward_sdsv;
    backward_sdsv.dye = 0.6;
    backward_sdsv.no_dye = 0.4;
    StuckDyeStateVector next_backward_sdsv;
    next_backward_sdsv.dye = 0.8;
    next_backward_sdsv.no_dye = 0.2;
    ErrorModelFitter emf;
    LogNormalDistributionFitter* original_dist_fit = emf.distribution_fit;
    Mock<LogNormalDistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    emf.distribution_fit = &df_mock.get();
    e.improve_fit(forward_sdsv,
                  backward_sdsv,
                  next_backward_sdsv,
                  num_edmans,
                  probability,
                  &emf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.2, TOL), 0, Close(0.3 * 0.4 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.2, TOL), 1, Close(0.7 * 0.6 / 0.98765, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    emf.distribution_fit = original_dist_fit;  // This avoids breaking cleanup.
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_timesteps_test, *tolerance(TOL)) {
    int num_timesteps = 2;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.2;
    rad(1, 0) = 1.3;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 1;
    double probability = 0.98765;
    StuckDyeStateVector forward_sdsv;
    forward_sdsv.dye = 0.7;
    forward_sdsv.no_dye = 0.3;
    StuckDyeStateVector backward_sdsv;
    backward_sdsv.dye = 0.6;
    backward_sdsv.no_dye = 0.4;
    StuckDyeStateVector next_backward_sdsv;
    next_backward_sdsv.dye = 0.8;
    next_backward_sdsv.no_dye = 0.2;
    ErrorModelFitter emf;
    LogNormalDistributionFitter* original_dist_fit = emf.distribution_fit;
    Mock<LogNormalDistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    emf.distribution_fit = &df_mock.get();
    e.improve_fit(forward_sdsv,
                  backward_sdsv,
                  next_backward_sdsv,
                  num_edmans,
                  probability,
                  &emf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.3, TOL), 0, Close(0.3 * 0.4 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.3, TOL), 1, Close(0.7 * 0.6 / 0.98765, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    emf.distribution_fit = original_dist_fit;  // This avoids breaking cleanup.
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_dye_colors_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.2;
    rad(0, 1) = 1.3;
    int channel = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.23 * (double)state;
    };
    StuckDyeEmission e(rad, channel, pdf);
    int num_edmans = 0;
    double probability = 0.98765;
    StuckDyeStateVector forward_sdsv;
    forward_sdsv.dye = 0.7;
    forward_sdsv.no_dye = 0.3;
    StuckDyeStateVector backward_sdsv;
    backward_sdsv.dye = 0.6;
    backward_sdsv.no_dye = 0.4;
    StuckDyeStateVector next_backward_sdsv;
    next_backward_sdsv.dye = 0.8;
    next_backward_sdsv.no_dye = 0.2;
    ErrorModelFitter emf;
    LogNormalDistributionFitter* original_dist_fit = emf.distribution_fit;
    Mock<LogNormalDistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    emf.distribution_fit = &df_mock.get();
    e.improve_fit(forward_sdsv,
                  backward_sdsv,
                  next_backward_sdsv,
                  num_edmans,
                  probability,
                  &emf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.2, TOL), 0, Close(0.3 * 0.4 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.2, TOL), 1, Close(0.7 * 0.6 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.3, TOL),
                          0,
                          Close((0.3 * 0.4 + 0.7 * 0.6) / 0.98765, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    emf.distribution_fit = original_dist_fit;  // This avoids breaking cleanup.
}

BOOST_AUTO_TEST_SUITE_END()  // stuck_dye_emission_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
