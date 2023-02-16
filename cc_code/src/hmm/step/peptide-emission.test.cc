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
#include "peptide-emission.h"

// Standard C++ library headers:
#include <functional>
#include <limits>

// External headers:
#include "fakeit.hpp"

// Local project headers:
#include "common/radiometry.h"
#include "parameterization/fit/log-normal-distribution-fitter.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "tensor/tensor.h"
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
BOOST_AUTO_TEST_SUITE(peptide_emission_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 2.0;
    int max_num_dyes = 5;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf)).AlwaysReturn(0.5);
    When(Method(cm_mock, sigma)).AlwaysReturn(0.6);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = 2.0;
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.num_channels == num_channels);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    BOOST_TEST((*e.values)[0] == 0.5);
    BOOST_TEST(e.pruned_range.min.size() == 2u);
    BOOST_TEST(e.pruned_range.min[0] == 0u);
    BOOST_TEST(e.pruned_range.min[1] == 1u);
    BOOST_TEST(e.pruned_range.max.size() == 2u);
    BOOST_TEST(e.pruned_range.max[0] == 1u);
    BOOST_TEST(e.pruned_range.max[1] == 4u);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(constructor_no_cutoff_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 2.0;
    int max_num_dyes = 5;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf)).AlwaysReturn(0.5);
    When(Method(cm_mock, sigma)).AlwaysReturn(0.6);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.num_channels == num_channels);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    BOOST_TEST((*e.values)[0] == 0.5);
    BOOST_TEST(e.pruned_range.min.size() == 2u);
    BOOST_TEST(e.pruned_range.min[0] == 0u);
    BOOST_TEST(e.pruned_range.min[1] == 0u);
    BOOST_TEST(e.pruned_range.max.size() == 2u);
    BOOST_TEST(e.pruned_range.max[0] == 1u);
    BOOST_TEST(e.pruned_range.max[1]
               == std::numeric_limits<unsigned int>::max());
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
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.prob(0, 0) == cm_mock.get().pdf(1.0, 0));
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
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    const PeptideEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0) == cm_mock.get().pdf(1.0, 0));
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
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.prob(0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(e.prob(1, 0) == cm_mock.get().pdf(0.1, 0));
    BOOST_TEST(e.prob(2, 0) == cm_mock.get().pdf(0.2, 0));
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
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    const PeptideEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(ce.prob(1, 0) == cm_mock.get().pdf(0.1, 0));
    BOOST_TEST(ce.prob(2, 0) == cm_mock.get().pdf(0.2, 0));
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
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock_0;
    When(Method(cm_mock_0, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock_0, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_0.get());
    Mock<ChannelModel> cm_mock_1;
    When(Method(cm_mock_1, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed * 1.9 + 0.098;
            });
    When(Method(cm_mock_1, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_1.get());
    Mock<ChannelModel> cm_mock_2;
    When(Method(cm_mock_2, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed * 2.7 + 0.187;
            });
    When(Method(cm_mock_2, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_2.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.prob(0, 0) == cm_mock_0.get().pdf(0.0, 0));
    BOOST_TEST(e.prob(1, 0) == cm_mock_1.get().pdf(0.1, 0));
    BOOST_TEST(e.prob(2, 0) == cm_mock_2.get().pdf(0.2, 0));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_dye_counts_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    int max_num_dyes = 2;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return 1.0 / (double)(state + 7);
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    BOOST_TEST(e.prob(0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 1) == cm_mock.get().pdf(0.0, 1));
    BOOST_TEST(e.prob(0, 2) == cm_mock.get().pdf(0.0, 2));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_dye_counts_const_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    int max_num_dyes = 2;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return 1.0 / (double)(state + 7);
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    const PeptideEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0) == cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 1) == cm_mock.get().pdf(0.0, 1));
    BOOST_TEST(ce.prob(0, 2) == cm_mock.get().pdf(0.0, 2));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_everything_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    int max_num_dyes = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return (observed + 0.042) / (double)(state + 7);
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    BOOST_TEST(e.prob(0, 0) == cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(e.prob(0, 1) == cm_mock.get().pdf(1.0, 1));
    BOOST_TEST(e.prob(1, 0) == cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(e.prob(1, 1) == cm_mock.get().pdf(1.1, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(prob_multiple_everything_const_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    int max_num_dyes = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return (observed + 0.042) / (double)(state + 7);
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int expected_size = num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values->size() == expected_size);
    const PeptideEmission& ce = e;
    BOOST_TEST(ce.prob(0, 0) == cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(ce.prob(0, 1) == cm_mock.get().pdf(1.0, 1));
    BOOST_TEST(ce.prob(1, 0) == cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(ce.prob(1, 1) == cm_mock.get().pdf(1.1, 1));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(forward_trivial_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return 0.5;
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 0};
    e.pruned_range.max = {num_timesteps, 1};
    e.allow_detached = true;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 3.14;
    psv1.broken_n_tensor[{0, 0}] = 3.014;
    psv1.allow_detached = true;
    psv1.p_detached = 1.23;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 3.14 * cm_mock.get().pdf(1.0, 0));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == 3.014 * cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(psv2->p_detached == 1.23 * cm_mock.get().pdf(1.0, 0));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 1u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_timesteps_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 2;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 0};
    e.pruned_range.max = {num_timesteps, 1};
    e.allow_detached = true;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 13.0;
    psv1.tensor[{1, 0}] = 13.1;
    psv1.tensor[{2, 0}] = 13.2;
    psv1.broken_n_tensor[{0, 0}] = 3.0;
    psv1.broken_n_tensor[{1, 0}] = 3.1;
    psv1.broken_n_tensor[{2, 0}] = 3.2;
    psv1.allow_detached = true;
    psv1.p_detached = 1.23;
    unsigned int edmans = 2;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 13.0 * cm_mock.get().pdf(2.0, 0));
    BOOST_TEST((psv2->tensor[{1, 0}]) == 13.1 * cm_mock.get().pdf(2.0, 0));
    BOOST_TEST((psv2->tensor[{2, 0}]) == 13.2 * cm_mock.get().pdf(2.0, 0));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == 3.0 * cm_mock.get().pdf(2.0, 0));
    BOOST_TEST((psv2->broken_n_tensor[{1, 0}])
               == 3.1 * cm_mock.get().pdf(2.0, 0));
    BOOST_TEST((psv2->broken_n_tensor[{2, 0}])
               == 3.2 * cm_mock.get().pdf(2.0, 0));
    BOOST_TEST(psv2->p_detached == 1.23 * cm_mock.get().pdf(2.0, 0));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.max[0] == 3u);
    BOOST_TEST(psv2->range.max[1] == 1u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 0, 0, 0};
    e.pruned_range.max = {num_timesteps, 1, 1, 1};
    e.allow_detached = true;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    shape[2] = 1;
    shape[3] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0, 0}] = 13.0;
    psv1.broken_n_tensor[{0, 0, 0, 0}] = 3.0;
    psv1.allow_detached = true;
    psv1.p_detached = 1.23;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0, 0, 0}])
               == 13.0 * cm_mock.get().pdf(0.0, 0) * cm_mock.get().pdf(0.1, 0)
                          * cm_mock.get().pdf(0.2, 0));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0, 0}])
               == 3.0 * cm_mock.get().pdf(0.0, 0) * cm_mock.get().pdf(0.1, 0)
                          * cm_mock.get().pdf(0.2, 0));
    BOOST_TEST(psv2->p_detached
               == 1.23 * cm_mock.get().pdf(0.0, 0) * cm_mock.get().pdf(0.1, 0)
                          * cm_mock.get().pdf(0.2, 0));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.min[2] == 0u);
    BOOST_TEST(psv2->range.min[3] == 0u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 1u);
    BOOST_TEST(psv2->range.max[2] == 1u);
    BOOST_TEST(psv2->range.max[3] == 1u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_different_pdfs_test,
                     *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock_0;
    When(Method(cm_mock_0, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed + 0.042;
            });
    When(Method(cm_mock_0, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_0.get());
    Mock<ChannelModel> cm_mock_1;
    When(Method(cm_mock_1, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed * 1.9 + 0.098;
            });
    When(Method(cm_mock_1, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_1.get());
    Mock<ChannelModel> cm_mock_2;
    When(Method(cm_mock_2, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return observed * 2.7 + 0.187;
            });
    When(Method(cm_mock_2, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_2.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 0, 0, 0};
    e.pruned_range.max = {num_timesteps, 1, 1, 1};
    e.allow_detached = true;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    shape[2] = 1;
    shape[3] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0, 0}] = 13.0;
    psv1.broken_n_tensor[{0, 0, 0, 0}] = 3.0;
    psv1.allow_detached = true;
    psv1.p_detached = 1.23;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0, 0, 0}])
               == 13.0 * cm_mock_0.get().pdf(0.0, 0)
                          * cm_mock_1.get().pdf(0.1, 0)
                          * cm_mock_2.get().pdf(0.2, 0));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0, 0}])
               == 3.0 * cm_mock_0.get().pdf(0.0, 0)
                          * cm_mock_1.get().pdf(0.1, 0)
                          * cm_mock_2.get().pdf(0.2, 0));
    BOOST_TEST(psv2->p_detached
               == 1.23 * cm_mock_0.get().pdf(0.0, 0)
                          * cm_mock_1.get().pdf(0.1, 0)
                          * cm_mock_2.get().pdf(0.2, 0));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.min[2] == 0u);
    BOOST_TEST(psv2->range.min[3] == 0u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 1u);
    BOOST_TEST(psv2->range.max[2] == 1u);
    BOOST_TEST(psv2->range.max[3] == 1u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_dye_counts_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    int max_num_dyes = 2;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return 1.0 / (double)(state + 7);
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 0};
    e.pruned_range.max = {num_timesteps, max_num_dyes + 1u};
    e.allow_detached = true;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 13.0;
    psv1.tensor[{0, 1}] = 13.1;
    psv1.tensor[{0, 2}] = 13.2;
    psv1.broken_n_tensor[{0, 0}] = 3.0;
    psv1.broken_n_tensor[{0, 1}] = 3.1;
    psv1.broken_n_tensor[{0, 2}] = 3.2;
    psv1.allow_detached = true;
    psv1.p_detached = 1.23;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == 13.0 * cm_mock.get().pdf(0.0, 0));
    BOOST_TEST((psv2->tensor[{0, 1}]) == 13.1 * cm_mock.get().pdf(0.0, 1));
    BOOST_TEST((psv2->tensor[{0, 2}]) == 13.2 * cm_mock.get().pdf(0.0, 2));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == 3.0 * cm_mock.get().pdf(0.0, 0));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == 3.1 * cm_mock.get().pdf(0.0, 1));
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}])
               == 3.2 * cm_mock.get().pdf(0.0, 2));
    BOOST_TEST(psv2->p_detached == 1.23 * cm_mock.get().pdf(0.0, 0));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 3u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_reduced_range_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    int max_num_dyes = 4;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return 1.0 / (double)(state + 7);
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 1};
    e.pruned_range.max = {num_timesteps, 3};
    e.allow_detached = false;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 1}] = 13.1;
    psv1.tensor[{0, 2}] = 13.2;
    psv1.broken_n_tensor[{0, 1}] = 3.1;
    psv1.broken_n_tensor[{0, 2}] = 3.2;
    psv1.allow_detached = false;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 1}]) == 13.1 * cm_mock.get().pdf(0.0, 1));
    BOOST_TEST((psv2->tensor[{0, 2}]) == 13.2 * cm_mock.get().pdf(0.0, 2));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == 3.1 * cm_mock.get().pdf(0.0, 1));
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}])
               == 3.2 * cm_mock.get().pdf(0.0, 2));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 1u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 3u);
    BOOST_TEST(psv2->allow_detached == false);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_multiple_everything_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    int max_num_dyes = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf))
            .AlwaysDo([](double observed, int state) -> double {
                return (observed + 0.042) / (double)(state + 7);
            });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 0, 0};
    e.pruned_range.max = {num_timesteps, 2, 2};
    e.allow_detached = true;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 7.000;
    psv1.tensor[{0, 0, 1}] = 7.001;
    psv1.tensor[{0, 1, 0}] = 7.010;
    psv1.tensor[{0, 1, 1}] = 7.011;
    psv1.tensor[{1, 0, 0}] = 7.100;
    psv1.tensor[{1, 0, 1}] = 7.101;
    psv1.tensor[{1, 1, 0}] = 7.110;
    psv1.tensor[{1, 1, 1}] = 7.111;
    psv1.broken_n_tensor[{0, 0, 0}] = 17.000;
    psv1.broken_n_tensor[{0, 0, 1}] = 17.001;
    psv1.broken_n_tensor[{0, 1, 0}] = 17.010;
    psv1.broken_n_tensor[{0, 1, 1}] = 17.011;
    psv1.broken_n_tensor[{1, 0, 0}] = 17.100;
    psv1.broken_n_tensor[{1, 0, 1}] = 17.101;
    psv1.broken_n_tensor[{1, 1, 0}] = 17.110;
    psv1.broken_n_tensor[{1, 1, 1}] = 17.111;
    psv1.allow_detached = true;
    psv1.p_detached = 1.23;
    unsigned int edmans = 1;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0, 0}])
               == 7.000 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->tensor[{0, 0, 1}])
               == 7.001 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST((psv2->tensor[{0, 1, 0}])
               == 7.010 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->tensor[{0, 1, 1}])
               == 7.011 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST((psv2->tensor[{1, 0, 0}])
               == 7.100 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->tensor[{1, 0, 1}])
               == 7.101 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST((psv2->tensor[{1, 1, 0}])
               == 7.110 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->tensor[{1, 1, 1}])
               == 7.111 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0}])
               == 17.000 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 1}])
               == 17.001 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 0}])
               == 17.010 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 1}])
               == 17.011 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST((psv2->broken_n_tensor[{1, 0, 0}])
               == 17.100 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->broken_n_tensor[{1, 0, 1}])
               == 17.101 * cm_mock.get().pdf(1.0, 0)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST((psv2->broken_n_tensor[{1, 1, 0}])
               == 17.110 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST((psv2->broken_n_tensor[{1, 1, 1}])
               == 17.111 * cm_mock.get().pdf(1.0, 1)
                          * cm_mock.get().pdf(1.1, 1));
    BOOST_TEST(psv2->p_detached
               == 1.23 * cm_mock.get().pdf(1.0, 0) * cm_mock.get().pdf(1.1, 0));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.min[2] == 0u);
    BOOST_TEST(psv2->range.max[0] == 2u);
    BOOST_TEST(psv2->range.max[1] == 2u);
    BOOST_TEST(psv2->range.max[2] == 2u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(improve_fit_simple_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.009;
    int max_num_dyes = 2;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf)).AlwaysReturn(0.5);
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 1.72;
    fpsv.tensor[{0, 1}] = 1.36;
    fpsv.tensor[{0, 2}] = 1.18;
    fpsv.broken_n_tensor[{0, 0}] = 1.072;
    fpsv.broken_n_tensor[{0, 1}] = 1.036;
    fpsv.broken_n_tensor[{0, 2}] = 1.018;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 2.72;
    bpsv.tensor[{0, 1}] = 2.36;
    bpsv.tensor[{0, 2}] = 2.18;
    bpsv.broken_n_tensor[{0, 0}] = 2.072;
    bpsv.broken_n_tensor[{0, 1}] = 2.036;
    bpsv.broken_n_tensor[{0, 2}] = 2.018;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 3.72;
    nbpsv.tensor[{0, 1}] = 3.36;
    nbpsv.tensor[{0, 2}] = 3.18;
    nbpsv.broken_n_tensor[{0, 0}] = 3.072;
    nbpsv.broken_n_tensor[{0, 1}] = 3.036;
    nbpsv.broken_n_tensor[{0, 2}] = 3.018;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 3.14159;
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    NormalDistributionFitter* original_dist_fit =
            smf.channel_fits[0]->distribution_fit;
    Mock<NormalDistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    smf.channel_fits[0]->distribution_fit = &df_mock.get();
    e.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.72 * 2.72 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.36 * 2.36 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          2,
                          Close(1.18 * 2.18 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.072 * 2.072 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.036 * 2.036 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          2,
                          Close(1.018 * 2.018 / 3.14159, TOL)))
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
    rad(0, 0) = 1.009;
    rad(0, 1) = 1.019;
    int max_num_dyes = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf)).AlwaysReturn(0.5);
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    shape[2] = max_num_dyes + 1;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0, 0}] = 1.72;
    fpsv.tensor[{0, 0, 1}] = 1.64;
    fpsv.tensor[{0, 1, 0}] = 1.36;
    fpsv.tensor[{0, 1, 1}] = 1.25;
    fpsv.broken_n_tensor[{0, 0, 0}] = 1.072;
    fpsv.broken_n_tensor[{0, 0, 1}] = 1.064;
    fpsv.broken_n_tensor[{0, 1, 0}] = 1.036;
    fpsv.broken_n_tensor[{0, 1, 1}] = 1.025;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0, 0}] = 2.72;
    bpsv.tensor[{0, 0, 1}] = 2.64;
    bpsv.tensor[{0, 1, 0}] = 2.36;
    bpsv.tensor[{0, 1, 1}] = 2.25;
    bpsv.broken_n_tensor[{0, 0, 0}] = 2.072;
    bpsv.broken_n_tensor[{0, 0, 1}] = 2.064;
    bpsv.broken_n_tensor[{0, 1, 0}] = 2.036;
    bpsv.broken_n_tensor[{0, 1, 1}] = 2.025;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0, 0}] = 3.72;
    nbpsv.tensor[{0, 0, 1}] = 3.64;
    nbpsv.tensor[{0, 1, 0}] = 3.36;
    nbpsv.tensor[{0, 1, 1}] = 3.25;
    nbpsv.broken_n_tensor[{0, 0, 0}] = 3.072;
    nbpsv.broken_n_tensor[{0, 0, 1}] = 3.064;
    nbpsv.broken_n_tensor[{0, 1, 0}] = 3.036;
    nbpsv.broken_n_tensor[{0, 1, 1}] = 3.025;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 3.14159;
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
    e.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.72 * 2.72 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          0,
                          Close(1.72 * 2.72 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.64 * 2.64 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          1,
                          Close(1.64 * 2.64 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.36 * 2.36 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          0,
                          Close(1.36 * 2.36 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.25 * 2.25 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          1,
                          Close(1.25 * 2.25 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.072 * 2.072 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          0,
                          Close(1.072 * 2.072 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.064 * 2.064 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          1,
                          Close(1.064 * 2.064 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.036 * 2.036 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          0,
                          Close(1.036 * 2.036 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_0, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.025 * 2.025 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock_1, add_sample)
                   .Using(Close(1.019, TOL),
                          1,
                          Close(1.025 * 2.025 / 3.14159, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock_0);
    VerifyNoOtherInvocations(df_mock_1);
    // This avoids breaking cleanup:
    smf.channel_fits[0]->distribution_fit = original_dist_fit_0;
    smf.channel_fits[1]->distribution_fit = original_dist_fit_1;
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_edmans_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.009;
    rad(1, 0) = 1.109;
    int max_num_dyes = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(Method(cm_mock, pdf)).AlwaysReturn(0.5);
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 1.72;
    fpsv.tensor[{0, 1}] = 1.36;
    fpsv.tensor[{1, 0}] = 1.64;
    fpsv.tensor[{1, 1}] = 1.25;
    fpsv.broken_n_tensor[{0, 0}] = 1.072;
    fpsv.broken_n_tensor[{0, 1}] = 1.036;
    fpsv.broken_n_tensor[{1, 0}] = 1.064;
    fpsv.broken_n_tensor[{1, 1}] = 1.025;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 2.72;
    bpsv.tensor[{0, 1}] = 2.36;
    bpsv.tensor[{1, 0}] = 2.64;
    bpsv.tensor[{1, 1}] = 2.25;
    bpsv.broken_n_tensor[{0, 0}] = 2.072;
    bpsv.broken_n_tensor[{0, 1}] = 2.036;
    bpsv.broken_n_tensor[{1, 0}] = 2.064;
    bpsv.broken_n_tensor[{1, 1}] = 2.025;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 3.72;
    nbpsv.tensor[{0, 1}] = 3.36;
    nbpsv.tensor[{1, 0}] = 3.64;
    nbpsv.tensor[{1, 1}] = 3.25;
    nbpsv.broken_n_tensor[{0, 0}] = 3.072;
    nbpsv.broken_n_tensor[{0, 1}] = 3.036;
    nbpsv.broken_n_tensor[{1, 0}] = 3.064;
    nbpsv.broken_n_tensor[{1, 1}] = 3.025;
    delete[] shape;
    unsigned int edmans = 1;
    double probability = 3.14159;
    SequencingModelFitter smf;
    smf.channel_fits.push_back(new ChannelModelFitter());
    NormalDistributionFitter* original_dist_fit =
            smf.channel_fits[0]->distribution_fit;
    Mock<NormalDistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    smf.channel_fits[0]->distribution_fit = &df_mock.get();
    e.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &smf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          0,
                          Close(1.72 * 2.72 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          1,
                          Close(1.36 * 2.36 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          0,
                          Close(1.64 * 2.64 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          1,
                          Close(1.25 * 2.25 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          0,
                          Close(1.072 * 2.072 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          1,
                          Close(1.036 * 2.036 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          0,
                          Close(1.064 * 2.064 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          1,
                          Close(1.025 * 2.025 / 3.14159, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    smf.channel_fits[0]->distribution_fit =
            original_dist_fit;  // This avoids breaking cleanup.
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_SUITE_END()  // peptide_emission_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
