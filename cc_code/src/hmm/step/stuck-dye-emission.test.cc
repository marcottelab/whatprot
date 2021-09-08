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
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "test-util/fakeit.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using fakeit::Fake;
using fakeit::Mock;
using fakeit::Verify;
using fakeit::VerifyNoOtherInvocations;
using fakeit::When;
using std::function;
using whatprot::test_util::Close;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(stuck_dye_emission_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.num_channels == num_channels);
    BOOST_TEST(e.channel == channel);
    BOOST_TEST(e.values[0] == 0.5);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_timesteps_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 0, 0) == cm_mock.get().pdf(0.0, 1));
    BOOST_TEST(e.prob(1, 0, 0) == cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(e.prob(1, 0, 0) == cm_mock.get().pdf(1.0, 1));
    BOOST_TEST(e.prob(2, 0, 0) == cm_mock.get().pdf(2.0, 0));
    BOOST_TEST(e.prob(2, 0, 0) == cm_mock.get().pdf(2.0, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_timesteps_const_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    const StuckDyeEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 0, 0) == cm_mock.get().pdf(0.0, 1));
    BOOST_TEST(ce.prob(1, 0, 0) == cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(ce.prob(1, 0, 0) == cm_mock.get().pdf(1.0, 1));
    BOOST_TEST(ce.prob(2, 0, 0) == cm_mock.get().pdf(2.0, 0));
    BOOST_TEST(ce.prob(2, 0, 0) == cm_mock.get().pdf(2.0, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_channels_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.27 * (double)state + 0.042;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 0, 1) == cm_mock.get().pdf(0.0, 1));
    BOOST_TEST(e.prob(0, 1, 0) == cm_mock.get().pdf(0.1, 0));
    BOOST_TEST(e.prob(0, 2, 0) == cm_mock.get().pdf(0.2, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_channels_different_pdfs_test,
                     *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock_0;
    When(Method(cm_mock_0, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.27 * (double)state + 0.042;
            });
    seq_model.channel_models.push_back(&cm_mock_0.get());
    Mock<ChannelModel> cm_mock_1;
    When(Method(cm_mock_1, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed * 1.9 + 0.34 * (double)state + 0.098;
            });
    seq_model.channel_models.push_back(&cm_mock_1.get());
    Mock<ChannelModel> cm_mock_2;
    When(Method(cm_mock_2, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed * 2.7 + 0.19 * (double)state + 0.187;
            });
    seq_model.channel_models.push_back(&cm_mock_2.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == cm_mock_0.get().pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 0, 1) == cm_mock_0.get().pdf(0.0, 1));
    BOOST_TEST(e.prob(0, 1, 0) == cm_mock_1.get().pdf(0.1, 0));
    BOOST_TEST(e.prob(0, 2, 0) == cm_mock_2.get().pdf(0.2, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_channels_const_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.27 * (double)state + 0.042;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int expected_size = num_timesteps * num_channels * 2;
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    const StuckDyeEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 0, 1) == cm_mock.get().pdf(0.0, 1));
    BOOST_TEST(ce.prob(0, 1, 0) == cm_mock.get().pdf(0.1, 0));
    BOOST_TEST(ce.prob(0, 2, 0) == cm_mock.get().pdf(0.2, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(forward_basic_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye == 0.3 * cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(sdsv.dye == 0.7 * cm_mock.get().pdf(1.0, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(forward_multiple_timestamps_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(1, 0) = 1.1;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 1;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye == 0.3 * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(sdsv.dye == 0.7 * cm_mock.get().pdf(1.1, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye
               == 0.3 * cm_mock.get().pdf(1.0, 0) * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(sdsv.dye
               == 0.7 * cm_mock.get().pdf(1.0, 1) * cm_mock.get().pdf(1.1, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_other_channel_test,
                     *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    unsigned int channel = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye
               == 0.3 * cm_mock.get().pdf(1.0, 0) * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(sdsv.dye
               == 0.7 * cm_mock.get().pdf(1.0, 0) * cm_mock.get().pdf(1.1, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_different_pdfs_test,
                     *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock_0;
    When(Method(cm_mock_0, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock_0.get());
    Mock<ChannelModel> cm_mock_1;
    When(Method(cm_mock_1, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return 2.1 * observed + 0.49 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock_1.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector sdsv;
    sdsv.dye = 0.7;
    sdsv.no_dye = 0.3;
    e.forward(&num_edmans, &sdsv);
    BOOST_TEST(sdsv.no_dye
               == 0.3 * cm_mock_0.get().pdf(1.0, 0)
                          * cm_mock_1.get().pdf(1.1, 0));
    BOOST_TEST(sdsv.dye
               == 0.7 * cm_mock_0.get().pdf(1.0, 1)
                          * cm_mock_1.get().pdf(1.1, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(backward_basic_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye == 0.3 * cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(output.dye == 0.7 * cm_mock.get().pdf(1.0, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(backward_multiple_timestamps_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(1, 0) = 1.1;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 1;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye == 0.3 * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(output.dye == 0.7 * cm_mock.get().pdf(1.1, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(backward_multiple_channels_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye
               == 0.3 * cm_mock.get().pdf(1.0, 0) * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(output.dye
               == 0.7 * cm_mock.get().pdf(1.0, 1) * cm_mock.get().pdf(1.1, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(backward_multiple_channels_other_channel_test,
                     *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    unsigned int channel = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye
               == 0.3 * cm_mock.get().pdf(1.0, 0) * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(output.dye
               == 0.7 * cm_mock.get().pdf(1.0, 0) * cm_mock.get().pdf(1.1, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(backward_multiple_channels_different_pdfs_test,
                     *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    rad(0, 1) = 1.1;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock_0;
    When(Method(cm_mock_0, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock_0.get());
    Mock<ChannelModel> cm_mock_1;
    When(Method(cm_mock_1, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return 2.1 * observed + 0.49 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock_1.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
    StuckDyeStateVector input;
    input.dye = 0.7;
    input.no_dye = 0.3;
    StuckDyeStateVector output;
    e.backward(input, &num_edmans, &output);
    BOOST_TEST(output.no_dye
               == 0.3 * cm_mock_0.get().pdf(1.0, 0)
                          * cm_mock_1.get().pdf(1.1, 0));
    BOOST_TEST(output.dye
               == 0.7 * cm_mock_0.get().pdf(1.0, 1)
                          * cm_mock_1.get().pdf(1.1, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(improve_fit_basic_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.2;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
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
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    NormalDistributionFitter* original_dist_fit =
            smf.channel_fits[0]->distribution_fit;
    Mock<NormalDistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    smf.channel_fits[0]->distribution_fit = &df_mock.get();
    e.improve_fit(forward_sdsv,
                  backward_sdsv,
                  next_backward_sdsv,
                  num_edmans,
                  probability,
                  &smf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.2, TOL), 0, Close(0.3 * 0.4 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.2, TOL), 1, Close(0.7 * 0.6 / 0.98765, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    smf.channel_fits[0]->distribution_fit =
            original_dist_fit;  // This avoids breaking cleanup.
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_timesteps_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.2;
    rad(1, 0) = 1.3;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 1;
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
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    NormalDistributionFitter* original_dist_fit =
            smf.channel_fits[0]->distribution_fit;
    Mock<NormalDistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    smf.channel_fits[0]->distribution_fit = &df_mock.get();
    e.improve_fit(forward_sdsv,
                  backward_sdsv,
                  next_backward_sdsv,
                  num_edmans,
                  probability,
                  &smf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.3, TOL), 0, Close(0.3 * 0.4 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.3, TOL), 1, Close(0.7 * 0.6 / 0.98765, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    smf.channel_fits[0]->distribution_fit =
            original_dist_fit;  // This avoids breaking cleanup.
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_dye_colors_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.2;
    rad(0, 1) = 1.3;
    unsigned int channel = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.23 * (double)state;
            });
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    StuckDyeEmission e(rad, channel, seq_model);
    unsigned int num_edmans = 0;
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
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    smf.channel_fits.push_back(new ChannelModelFitter());
    NormalDistributionFitter* original_dist_fit_0 =
            smf.channel_fits[0]->distribution_fit;
    Mock<NormalDistributionFitter> df_mock_0;
    Fake(Method(df_mock_0, add_sample));
    smf.channel_fits[0]->distribution_fit = &df_mock_0.get();
    NormalDistributionFitter* original_dist_fit_1 =
            smf.channel_fits[1]->distribution_fit;
    Mock<NormalDistributionFitter> df_mock_1;
    Fake(Method(df_mock_1, add_sample));
    smf.channel_fits[1]->distribution_fit = &df_mock_1.get();
    e.improve_fit(forward_sdsv,
                  backward_sdsv,
                  next_backward_sdsv,
                  num_edmans,
                  probability,
                  &smf);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.2, TOL), 0, Close(0.3 * 0.4 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.2, TOL), 1, Close(0.7 * 0.6 / 0.98765, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.3, TOL),
                          0,
                          Close((0.3 * 0.4 + 0.7 * 0.6) / 0.98765, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock_0);
    VerifyNoOtherInvocations(df_mock_1);
    // This avoids breaking cleanup:
    smf.channel_fits[0]->distribution_fit = original_dist_fit_0;
    smf.channel_fits[1]->distribution_fit = original_dist_fit_1;
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_SUITE_END()  // stuck_dye_emission_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
