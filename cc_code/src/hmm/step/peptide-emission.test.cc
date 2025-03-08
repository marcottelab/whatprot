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
#include <vector>

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
using std::numeric_limits;
using std::vector;
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
    unsigned int max_num_dyes = 5;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysReturn(0.5);
    When(Method(cm_mock, sigma)).AlwaysReturn(0.6);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    BOOST_TEST(e.ptsr->range.min.size() == num_channels);
    BOOST_TEST(e.ptsr->range.min[0] == 0u);
    BOOST_TEST(e.ptsr->range.max.size() == num_channels);
    BOOST_TEST(e.ptsr->range.max[0] == max_num_dyes + 1);
    BOOST_TEST(e.num_channels == num_channels);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    BOOST_TEST(e.ptsr->values[0] == 0.5);
    BOOST_TEST(e.pruned_range.min.size() == 2u);
    BOOST_TEST(e.pruned_range.min[0] == 0u);
    BOOST_TEST(e.pruned_range.min[1] == 0u);
    BOOST_TEST(e.pruned_range.max.size() == 2u);
    BOOST_TEST(e.pruned_range.max[0] == 1u);
    BOOST_TEST(e.pruned_range.max[1] == max_num_dyes + 1);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(constructor_no_cutoff_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 2.0;
    unsigned int max_num_dyes = 5;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysReturn(0.5);
    When(Method(cm_mock, sigma)).AlwaysReturn(0.6);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    BOOST_TEST(e.ptsr->range.min.size() == num_channels);
    BOOST_TEST(e.ptsr->range.min[0] == 0u);
    BOOST_TEST(e.ptsr->range.max.size() == num_channels);
    BOOST_TEST(e.ptsr->range.max[0] == max_num_dyes + 1);
    BOOST_TEST(e.num_channels == num_channels);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    BOOST_TEST(e.ptsr->values[0] == 0.5);
    BOOST_TEST(e.pruned_range.min.size() == 2u);
    BOOST_TEST(e.pruned_range.min[0] == 0u);
    BOOST_TEST(e.pruned_range.min[1] == 0u);
    BOOST_TEST(e.pruned_range.max.size() == 2u);
    BOOST_TEST(e.pruned_range.max[0] == 1u);
    BOOST_TEST(e.pruned_range.max[1] == max_num_dyes + 1);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(ptsr_multiple_timesteps_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return observed + 0.042;
                    });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    BOOST_TEST(e.ptsr->range.min.size() == num_channels);
    BOOST_TEST(e.ptsr->range.min[0] == 0u);
    BOOST_TEST(e.ptsr->range.max.size() == num_channels);
    BOOST_TEST(e.ptsr->range.max[0] == max_num_dyes + 1);
    vector<unsigned int> counts(num_channels);
    counts = {0};
    BOOST_TEST((*e.ptsr)[&counts[0]] == cm_mock.get().pdf(1.0, &counts[0]));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(ptsr_multiple_channels_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
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
    BOOST_TEST(e.ptsr->range.min.size() == num_channels);
    BOOST_TEST(e.ptsr->range.min[0] == 0u);
    BOOST_TEST(e.ptsr->range.min[1] == 0u);
    BOOST_TEST(e.ptsr->range.min[2] == 0u);
    BOOST_TEST(e.ptsr->range.max.size() == num_channels);
    BOOST_TEST(e.ptsr->range.max[0] == max_num_dyes + 1);
    BOOST_TEST(e.ptsr->range.max[1] == max_num_dyes + 1);
    BOOST_TEST(e.ptsr->range.max[2] == max_num_dyes + 1);
    vector<unsigned int> counts(num_channels);
    counts = {0, 0, 0};
    BOOST_TEST((*e.ptsr)[&counts[0]]
               == cm_mock.get().pdf(0.0, &counts[0])
                          * cm_mock.get().pdf(0.1, &counts[0])
                          * cm_mock.get().pdf(0.2, &counts[0]));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(ptsr_multiple_channels_different_pdfs_test,
                     *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock_0;
    When(ConstOverloadedMethod(
                 cm_mock_0, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return observed + 0.042;
                    });
    When(Method(cm_mock_0, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_0.get());
    Mock<ChannelModel> cm_mock_1;
    When(ConstOverloadedMethod(
                 cm_mock_1, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return observed * 1.9 + 0.098;
                    });
    When(Method(cm_mock_1, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_1.get());
    Mock<ChannelModel> cm_mock_2;
    When(ConstOverloadedMethod(
                 cm_mock_2, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return observed * 2.7 + 0.187;
                    });
    When(Method(cm_mock_2, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_2.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    BOOST_TEST(e.ptsr->range.min.size() == num_channels);
    BOOST_TEST(e.ptsr->range.min[0] == 0u);
    BOOST_TEST(e.ptsr->range.min[1] == 0u);
    BOOST_TEST(e.ptsr->range.min[2] == 0u);
    BOOST_TEST(e.ptsr->range.max.size() == num_channels);
    BOOST_TEST(e.ptsr->range.max[0] == max_num_dyes + 1);
    BOOST_TEST(e.ptsr->range.max[1] == max_num_dyes + 1);
    BOOST_TEST(e.ptsr->range.max[2] == max_num_dyes + 1);
    vector<unsigned int> counts(num_channels);
    counts = {0, 0, 0};
    BOOST_TEST((*e.ptsr)[&counts[0]]
               == cm_mock_0.get().pdf(0.0, &counts[0])
                          * cm_mock_1.get().pdf(0.1, &counts[0])
                          * cm_mock_2.get().pdf(0.2, &counts[0]));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(ptsr_multiple_dye_counts_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    unsigned int max_num_dyes = 2;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return 1.0 / (double)(counts[0] + 7);
                    });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    BOOST_TEST(e.ptsr->range.min.size() == num_channels);
    BOOST_TEST(e.ptsr->range.min[0] == 0u);
    BOOST_TEST(e.ptsr->range.max.size() == num_channels);
    BOOST_TEST(e.ptsr->range.max[0] == max_num_dyes + 1);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    vector<unsigned int> counts(num_channels);
    counts = {0};
    BOOST_TEST((*e.ptsr)[&counts[0]] == cm_mock.get().pdf(0.0, &counts[0]));
    counts = {1};
    BOOST_TEST((*e.ptsr)[&counts[0]] == cm_mock.get().pdf(0.0, &counts[0]));
    counts = {2};
    BOOST_TEST((*e.ptsr)[&counts[0]] == cm_mock.get().pdf(0.0, &counts[0]));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(ptsr_multiple_everything_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 2;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    unsigned int max_num_dyes = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return (observed + 0.042)
                               / (double)(counts[0] + 3.14 * counts[1] + 7);
                    });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 1;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    BOOST_TEST(e.ptsr->range.min.size() == num_channels);
    BOOST_TEST(e.ptsr->range.min[0] == 0u);
    BOOST_TEST(e.ptsr->range.min[1] == 0u);
    BOOST_TEST(e.ptsr->range.max.size() == num_channels);
    BOOST_TEST(e.ptsr->range.max[0] == max_num_dyes + 1);
    BOOST_TEST(e.ptsr->range.max[1] == max_num_dyes + 1);
    vector<unsigned int> counts(num_channels);
    counts = {0, 0};
    BOOST_TEST((*e.ptsr)[&counts[0]]
               == cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 1};
    BOOST_TEST((*e.ptsr)[&counts[0]]
               == cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 0};
    BOOST_TEST((*e.ptsr)[&counts[0]]
               == cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 1};
    BOOST_TEST((*e.ptsr)[&counts[0]]
               == cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
}

BOOST_AUTO_TEST_CASE(forward_trivial_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
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
    vector<unsigned int> counts(num_channels);
    counts = {0};
    BOOST_TEST((psv2->tensor[{0, 0}])
               == 3.14 * cm_mock.get().pdf(1.0, &counts[0]));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == 3.014 * cm_mock.get().pdf(1.0, &counts[0]));
    BOOST_TEST(psv2->p_detached == 1.23 * cm_mock.get().pdf(1.0, &counts[0]));
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
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
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
    vector<unsigned int> counts(num_channels);
    counts = {0};
    BOOST_TEST((psv2->tensor[{0, 0}])
               == 13.0 * cm_mock.get().pdf(2.0, &counts[0]));
    BOOST_TEST((psv2->tensor[{1, 0}])
               == 13.1 * cm_mock.get().pdf(2.0, &counts[0]));
    BOOST_TEST((psv2->tensor[{2, 0}])
               == 13.2 * cm_mock.get().pdf(2.0, &counts[0]));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == 3.0 * cm_mock.get().pdf(2.0, &counts[0]));
    BOOST_TEST((psv2->broken_n_tensor[{1, 0}])
               == 3.1 * cm_mock.get().pdf(2.0, &counts[0]));
    BOOST_TEST((psv2->broken_n_tensor[{2, 0}])
               == 3.2 * cm_mock.get().pdf(2.0, &counts[0]));
    BOOST_TEST(psv2->p_detached == 1.23 * cm_mock.get().pdf(2.0, &counts[0]));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.max[0] == 3u);
    BOOST_TEST(psv2->range.max[1] == 1u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_two_channels_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return (observed + 0.042) / (counts[0] * 3 + counts[1]);
                    });
    When(Method(cm_mock, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock.get());
    seq_model.channel_models.push_back(&cm_mock.get());
    unsigned int timestep = 0;
    SequencingSettings seq_settings;
    seq_settings.dist_cutoff = std::numeric_limits<double>::max();
    PeptideEmission e(rad, timestep, max_num_dyes, seq_model, seq_settings);
    e.pruned_range.min = {0, 0, 0};
    e.pruned_range.max = {num_timesteps, 1, 1};
    e.allow_detached = true;
    unsigned int order = 1 + num_channels;
    unsigned int* shape = new unsigned int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    shape[2] = 1;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0}] = 13.0;
    psv1.broken_n_tensor[{0, 0, 0}] = 3.0;
    psv1.allow_detached = true;
    psv1.p_detached = 1.23;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = e.forward(psv1, &edmans);
    vector<unsigned int> counts(num_channels);
    counts = {0, 0};
    BOOST_TEST((psv2->tensor[{0, 0, 0}])
               == 13.0 * cm_mock.get().pdf(0.0, &counts[0])
                          * cm_mock.get().pdf(0.1, &counts[0]));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0}])
               == 3.0 * cm_mock.get().pdf(0.0, &counts[0])
                          * cm_mock.get().pdf(0.1, &counts[0]));
    BOOST_TEST(psv2->p_detached
               == 1.23 * cm_mock.get().pdf(0.0, &counts[0])
                          * cm_mock.get().pdf(0.1, &counts[0]));
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.min[2] == 0u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 1u);
    BOOST_TEST(psv2->range.max[2] == 1u);
    BOOST_TEST(psv2->allow_detached == true);
    // Avoid double clean-up:
    seq_model.channel_models.resize(0);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(forward_three_channels_test, *tolerance(TOL)) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
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
    vector<unsigned int> counts(num_channels);
    counts = {0, 0, 0};
    BOOST_TEST((psv2->tensor[{0, 0, 0, 0}])
               == 13.0 * cm_mock.get().pdf(0.0, &counts[0])
                          * cm_mock.get().pdf(0.1, &counts[0])
                          * cm_mock.get().pdf(0.2, &counts[0]));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0, 0}])
               == 3.0 * cm_mock.get().pdf(0.0, &counts[0])
                          * cm_mock.get().pdf(0.1, &counts[0])
                          * cm_mock.get().pdf(0.2, &counts[0]));
    BOOST_TEST(psv2->p_detached
               == 1.23 * cm_mock.get().pdf(0.0, &counts[0])
                          * cm_mock.get().pdf(0.1, &counts[0])
                          * cm_mock.get().pdf(0.2, &counts[0]));
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
    unsigned int max_num_dyes = 0;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock_0;
    When(ConstOverloadedMethod(
                 cm_mock_0, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return observed + 0.042;
                    });
    When(Method(cm_mock_0, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_0.get());
    Mock<ChannelModel> cm_mock_1;
    When(ConstOverloadedMethod(
                 cm_mock_1, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return observed * 1.9 + 0.098;
                    });
    When(Method(cm_mock_1, sigma)).AlwaysReturn(0.5);
    seq_model.channel_models.push_back(&cm_mock_1.get());
    Mock<ChannelModel> cm_mock_2;
    When(ConstOverloadedMethod(
                 cm_mock_2, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
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
    vector<unsigned int> counts(num_channels);
    counts = {0, 0, 0};
    BOOST_TEST((psv2->tensor[{0, 0, 0, 0}])
               == 13.0 * cm_mock_0.get().pdf(0.0, &counts[0])
                          * cm_mock_1.get().pdf(0.1, &counts[0])
                          * cm_mock_2.get().pdf(0.2, &counts[0]));
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0, 0}])
               == 3.0 * cm_mock_0.get().pdf(0.0, &counts[0])
                          * cm_mock_1.get().pdf(0.1, &counts[0])
                          * cm_mock_2.get().pdf(0.2, &counts[0]));
    BOOST_TEST(psv2->p_detached
               == 1.23 * cm_mock_0.get().pdf(0.0, &counts[0])
                          * cm_mock_1.get().pdf(0.1, &counts[0])
                          * cm_mock_2.get().pdf(0.2, &counts[0]));
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
    unsigned int max_num_dyes = 2;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return 1.0 / (double)(counts[0] + 7);
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
    vector<unsigned int> counts(num_channels);
    counts = {0};
    BOOST_TEST((psv2->tensor[{0, 0}])
               == 13.0 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {1};
    BOOST_TEST((psv2->tensor[{0, 1}])
               == 13.1 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {2};
    BOOST_TEST((psv2->tensor[{0, 2}])
               == 13.2 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {0};
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}])
               == 3.0 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {1};
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == 3.1 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {2};
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}])
               == 3.2 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {0};
    BOOST_TEST(psv2->p_detached == 1.23 * cm_mock.get().pdf(0.0, &counts[0]));
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
    unsigned int max_num_dyes = 4;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return 1.0 / (double)(counts[0] + 7);
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
    vector<unsigned int> counts(num_channels);
    counts = {1};
    BOOST_TEST((psv2->tensor[{0, 1}])
               == 13.1 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {2};
    BOOST_TEST((psv2->tensor[{0, 2}])
               == 13.2 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {1};
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}])
               == 3.1 * cm_mock.get().pdf(0.0, &counts[0]));
    counts = {2};
    BOOST_TEST((psv2->broken_n_tensor[{0, 2}])
               == 3.2 * cm_mock.get().pdf(0.0, &counts[0]));
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
    unsigned int max_num_dyes = 1;
    SequencingModel seq_model;
    Mock<ChannelModel> cm_mock;
    When(ConstOverloadedMethod(
                 cm_mock, pdf, double(double, const unsigned int*)))
            .AlwaysDo(
                    [](double observed, const unsigned int* counts) -> double {
                        return (observed + 0.042)
                               / (double)(counts[0] + 3 * counts[1] + 7);
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
    vector<unsigned int> counts(num_channels);
    counts = {0, 0};
    BOOST_TEST((psv2->tensor[{0, 0, 0}])
               == 7.000 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 1};
    BOOST_TEST((psv2->tensor[{0, 0, 1}])
               == 7.001 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 0};
    BOOST_TEST((psv2->tensor[{0, 1, 0}])
               == 7.010 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 1};
    BOOST_TEST((psv2->tensor[{0, 1, 1}])
               == 7.011 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 0};
    BOOST_TEST((psv2->tensor[{1, 0, 0}])
               == 7.100 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 1};
    BOOST_TEST((psv2->tensor[{1, 0, 1}])
               == 7.101 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 0};
    BOOST_TEST((psv2->tensor[{1, 1, 0}])
               == 7.110 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 1};
    BOOST_TEST((psv2->tensor[{1, 1, 1}])
               == 7.111 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 0};
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 0}])
               == 17.000 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 1};
    BOOST_TEST((psv2->broken_n_tensor[{0, 0, 1}])
               == 17.001 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 0};
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 0}])
               == 17.010 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 1};
    BOOST_TEST((psv2->broken_n_tensor[{0, 1, 1}])
               == 17.011 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 0};
    BOOST_TEST((psv2->broken_n_tensor[{1, 0, 0}])
               == 17.100 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 1};
    BOOST_TEST((psv2->broken_n_tensor[{1, 0, 1}])
               == 17.101 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 0};
    BOOST_TEST((psv2->broken_n_tensor[{1, 1, 0}])
               == 17.110 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {1, 1};
    BOOST_TEST((psv2->broken_n_tensor[{1, 1, 1}])
               == 17.111 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
    counts = {0, 0};
    BOOST_TEST(psv2->p_detached
               == 1.23 * cm_mock.get().pdf(1.0, &counts[0])
                          * cm_mock.get().pdf(1.1, &counts[0]));
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

BOOST_AUTO_TEST_SUITE_END()  // peptide_emission_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
