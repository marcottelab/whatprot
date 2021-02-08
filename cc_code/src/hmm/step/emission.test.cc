/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// These tests may benefit from using a mocking framework. If a mocking
// framework is added to this code base in the future many of these tests should
// be revisited with that in mind.

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// File under test:
#include "emission.h"

// Standard C++ library headers:
#include <functional>

// External headers:
#include "fakeit.hpp"

// Local project headers:
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/fit/distribution-fitter.h"
#include "hmm/fit/error-model-fitter.h"
#include "tensor/tensor.h"
#include "test-util/fakeit.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
using fakeit::Fake;
using fakeit::Mock;
using fakeit::Verify;
using fakeit::VerifyNoOtherInvocations;
using whatprot::test_util::Close;
using std::function;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(emission_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 0.5;
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.num_channels == num_channels);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    BOOST_TEST(e.values[0] == 0.5);
}

BOOST_AUTO_TEST_CASE(multiple_timesteps_test, *tolerance(TOL)) {
    int num_timesteps = 3;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(e.prob(1, 0, 0) == pdf(1.0, 0));
    BOOST_TEST(e.prob(2, 0, 0) == pdf(2.0, 0));
}

BOOST_AUTO_TEST_CASE(multiple_timesteps_const_test, *tolerance(TOL)) {
    int num_timesteps = 3;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    const Emission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(ce.prob(1, 0, 0) == pdf(1.0, 0));
    BOOST_TEST(ce.prob(2, 0, 0) == pdf(2.0, 0));
}

BOOST_AUTO_TEST_CASE(multiple_channels_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 1, 0) == pdf(0.1, 0));
    BOOST_TEST(e.prob(0, 2, 0) == pdf(0.2, 0));
}

BOOST_AUTO_TEST_CASE(multiple_channels_const_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    const Emission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 1, 0) == pdf(0.1, 0));
    BOOST_TEST(ce.prob(0, 2, 0) == pdf(0.2, 0));
}

BOOST_AUTO_TEST_CASE(multiple_dye_counts_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    int max_num_dyes = 2;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0 / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 0, 1) == pdf(0.0, 1));
    BOOST_TEST(e.prob(0, 0, 2) == pdf(0.0, 2));
}

BOOST_AUTO_TEST_CASE(multiple_dye_counts_const_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    int max_num_dyes = 2;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0 / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.max_num_dyes == max_num_dyes);
    const Emission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 0, 1) == pdf(0.0, 1));
    BOOST_TEST(ce.prob(0, 0, 2) == pdf(0.0, 2));
}

BOOST_AUTO_TEST_CASE(multiple_everything_test, *tolerance(TOL)) {
    int num_timesteps = 2;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    int max_num_dyes = 1;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return (observed + 0.042) / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    BOOST_TEST(e.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(e.prob(0, 0, 1) == pdf(0.0, 1));
    BOOST_TEST(e.prob(0, 1, 0) == pdf(0.1, 0));
    BOOST_TEST(e.prob(0, 1, 1) == pdf(0.1, 1));
    BOOST_TEST(e.prob(1, 0, 0) == pdf(1.0, 0));
    BOOST_TEST(e.prob(1, 0, 1) == pdf(1.0, 1));
    BOOST_TEST(e.prob(1, 1, 0) == pdf(1.1, 0));
    BOOST_TEST(e.prob(1, 1, 1) == pdf(1.1, 1));
}

BOOST_AUTO_TEST_CASE(multiple_everything_const_test, *tolerance(TOL)) {
    int num_timesteps = 2;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    int max_num_dyes = 1;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return (observed + 0.042) / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int expected_size = num_timesteps * num_channels * (max_num_dyes + 1);
    BOOST_TEST(e.values.size() == expected_size);
    BOOST_TEST(e.num_timesteps == num_timesteps);
    const Emission& ce = e;
    BOOST_TEST(ce.prob(0, 0, 0) == pdf(0.0, 0));
    BOOST_TEST(ce.prob(0, 0, 1) == pdf(0.0, 1));
    BOOST_TEST(ce.prob(0, 1, 0) == pdf(0.1, 0));
    BOOST_TEST(ce.prob(0, 1, 1) == pdf(0.1, 1));
    BOOST_TEST(ce.prob(1, 0, 0) == pdf(1.0, 0));
    BOOST_TEST(ce.prob(1, 0, 1) == pdf(1.0, 1));
    BOOST_TEST(ce.prob(1, 1, 0) == pdf(1.1, 0));
    BOOST_TEST(ce.prob(1, 1, 1) == pdf(1.1, 1));
}

BOOST_AUTO_TEST_CASE(forward_trivial_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 0.5;
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 3.14;  // loc is {0, 0}
    int edmans = 0;
    e.forward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 3.14 * pdf(1.0, 0));  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_tensor_reuse_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.0;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 0.5;
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 3.14;  // loc is {0, 0}
    int edmans = 0;
    e.forward(tsr, &edmans, &tsr);
    e.forward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 3.14 * pdf(1.0, 0) * pdf(1.0, 0));  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_multiple_timesteps_test, *tolerance(TOL)) {
    int num_timesteps = 3;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(1, 0) = 1.0;
    rad(2, 0) = 2.0;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 13.0;  // loc is {0, 0}
    loc[0] = 1;
    tsr[loc] = 13.1;  // loc is {1, 0}
    loc[0] = 2;
    tsr[loc] = 13.2;  // loc is {2, 0}
    int edmans = 2;
    e.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    BOOST_TEST(tsr[loc] == 13.0 * pdf(2.0, 0));  // loc is {0, 0}
    loc[0] = 1;
    BOOST_TEST(tsr[loc] == 13.1 * pdf(2.0, 0));  // loc is {1, 0}
    loc[0] = 2;
    BOOST_TEST(tsr[loc] == 13.2 * pdf(2.0, 0));  // loc is {2, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_multiple_channels_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 3;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(0, 2) = 0.2;
    int max_num_dyes = 0;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return observed + 0.042;
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = 1;
    shape[2] = 1;
    shape[3] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    tsr[loc] = 13.0;  // loc is {0, 0, 0, 0}
    int edmans = 0;
    e.forward(tsr, &edmans, &tsr);
    // loc is {0, 0, 0, 0}
    BOOST_TEST(tsr[loc] == 13.0 * pdf(0.0, 0) * pdf(0.1, 0) * pdf(0.2, 0));
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_multiple_dye_counts_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    int max_num_dyes = 2;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0 / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 13.0;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 13.1;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 13.2;  // loc is {0, 2}
    int edmans = 0;
    e.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 13.0 * pdf(0.0, 0));  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 13.1 * pdf(0.0, 1));  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 13.2 * pdf(0.0, 2));  // loc is {0, 2}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_multiple_everything_test, *tolerance(TOL)) {
    int num_timesteps = 2;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 0.0;
    rad(0, 1) = 0.1;
    rad(1, 0) = 1.0;
    rad(1, 1) = 1.1;
    int max_num_dyes = 1;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return (observed + 0.042) / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 7.000;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = 7.001;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 7.010;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = 7.011;  // loc is {0, 1, 1}
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 7.100;  // loc is {1, 0, 0}
    loc[2] = 1;
    tsr[loc] = 7.101;  // loc is {1, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 7.110;  // loc is {1, 1, 0}
    loc[2] = 1;
    tsr[loc] = 7.111;  // loc is {1, 1, 1}
    int edmans = 1;
    e.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {0, 0, 0}
    BOOST_TEST(tsr[loc] == 7.000 * pdf(1.0, 0) * pdf(1.1, 0));
    loc[2] = 1;
    // loc is {0, 0, 1}
    BOOST_TEST(tsr[loc] == 7.001 * pdf(1.0, 0) * pdf(1.1, 1));
    loc[1] = 1;
    loc[2] = 0;
    // loc is {0, 1, 0}
    BOOST_TEST(tsr[loc] == 7.010 * pdf(1.0, 1) * pdf(1.1, 0));
    loc[2] = 1;
    // loc is {0, 1, 1}
    BOOST_TEST(tsr[loc] == 7.011 * pdf(1.0, 1) * pdf(1.1, 1));
    loc[0] = 1;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {1, 0, 0}
    BOOST_TEST(tsr[loc] == 7.100 * pdf(1.0, 0) * pdf(1.1, 0));
    loc[2] = 1;
    // loc is {1, 0, 1}
    BOOST_TEST(tsr[loc] == 7.101 * pdf(1.0, 0) * pdf(1.1, 1));
    loc[1] = 1;
    loc[2] = 0;
    // loc is {1, 1, 0}
    BOOST_TEST(tsr[loc] == 7.110 * pdf(1.0, 1) * pdf(1.1, 0));
    loc[2] = 1;
    // loc is {1, 1, 1}
    BOOST_TEST(tsr[loc] == 7.111 * pdf(1.0, 1) * pdf(1.1, 1));
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(improve_fit_simple_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.009;
    int max_num_dyes = 2;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0 / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    Tensor ftsr(order, shape);
    Tensor btsr(order, shape);
    Tensor nbtsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    ftsr[loc] = 1.72;  // loc is {0, 0}
    btsr[loc] = 2.72;
    nbtsr[loc] = 3.72;
    loc[1] = 1;
    ftsr[loc] = 1.36;  // loc is {0, 1}
    btsr[loc] = 2.36;
    nbtsr[loc] = 3.36;
    loc[1] = 2;
    ftsr[loc] = 1.18;  // loc is {0, 2}
    btsr[loc] = 2.18;
    nbtsr[loc] = 3.18;
    delete[] loc;
    int edmans = 0;
    double probability = 3.14159;
    ErrorModelFitter emf(DistributionType::LOGNORMAL);
    DistributionFitter* original_dist_fit = emf.distribution_fit;
    Mock<DistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    emf.distribution_fit = &df_mock.get();
    e.improve_fit(ftsr, btsr, nbtsr, edmans, probability, &emf);
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
    VerifyNoOtherInvocations(df_mock);
    emf.distribution_fit = original_dist_fit;  // This avoids breaking cleanup.
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_dye_colors_test, *tolerance(TOL)) {
    int num_timesteps = 1;
    int num_channels = 2;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.009;
    rad(0, 1) = 1.019;
    int max_num_dyes = 1;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0 / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    shape[2] = max_num_dyes + 1;
    Tensor ftsr(order, shape);
    Tensor btsr(order, shape);
    Tensor nbtsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    ftsr[loc] = 1.72;  // loc is {0, 0, 0}
    btsr[loc] = 2.72;
    nbtsr[loc] = 3.72;
    loc[2] = 1;
    ftsr[loc] = 1.64;  // loc is {0, 0, 1}
    btsr[loc] = 2.64;
    nbtsr[loc] = 3.64;
    loc[1] = 1;
    loc[2] = 0;
    ftsr[loc] = 1.36;  // loc is {0, 1, 0}
    btsr[loc] = 2.36;
    nbtsr[loc] = 3.36;
    loc[2] = 1;
    ftsr[loc] = 1.25;  // loc is {0, 1, 1}
    btsr[loc] = 2.25;
    nbtsr[loc] = 3.25;
    delete[] loc;
    int edmans = 0;
    double probability = 3.14159;
    ErrorModelFitter emf(DistributionType::LOGNORMAL);
    DistributionFitter* original_dist_fit = emf.distribution_fit;
    Mock<DistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    emf.distribution_fit = &df_mock.get();
    e.improve_fit(ftsr, btsr, nbtsr, edmans, probability, &emf);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.72 * 2.72 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.019, TOL),
                          0,
                          Close(1.72 * 2.72 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          0,
                          Close(1.64 * 2.64 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.019, TOL),
                          1,
                          Close(1.64 * 2.64 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.36 * 2.36 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.019, TOL),
                          0,
                          Close(1.36 * 2.36 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.009, TOL),
                          1,
                          Close(1.25 * 2.25 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.019, TOL),
                          1,
                          Close(1.25 * 2.25 / 3.14159, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    emf.distribution_fit = original_dist_fit;  // This avoids breaking cleanup.
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_edmans_test, *tolerance(TOL)) {
    int num_timesteps = 2;
    int num_channels = 1;
    Radiometry rad(num_timesteps, num_channels);
    rad(0, 0) = 1.009;
    rad(1, 0) = 1.109;
    int max_num_dyes = 1;
    function<double(double, int)> pdf = [](double observed,
                                           int state) -> double {
        return 1.0 / (double)(state + 7);
    };
    Emission e(rad, max_num_dyes, pdf);
    int order = 1 + num_channels;
    int* shape = new int[order];
    shape[0] = num_timesteps;
    shape[1] = max_num_dyes + 1;
    Tensor ftsr(order, shape);
    Tensor btsr(order, shape);
    Tensor nbtsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    ftsr[loc] = 1.72;  // loc is {0, 0}
    btsr[loc] = 2.72;
    nbtsr[loc] = 3.72;
    loc[1] = 1;
    ftsr[loc] = 1.36;  // loc is {0, 1}
    btsr[loc] = 2.36;
    nbtsr[loc] = 3.36;
    loc[0] = 1;
    loc[1] = 0;
    ftsr[loc] = 1.64;  // loc is {1, 0}
    btsr[loc] = 2.64;
    nbtsr[loc] = 3.64;
    loc[1] = 1;
    ftsr[loc] = 1.25;  // loc is {1, 1}
    btsr[loc] = 2.25;
    nbtsr[loc] = 3.25;
    delete[] loc;
    int edmans = 1;
    double probability = 3.14159;
    ErrorModelFitter emf(DistributionType::LOGNORMAL);
    DistributionFitter* original_dist_fit = emf.distribution_fit;
    Mock<DistributionFitter> df_mock;
    Fake(Method(df_mock, add_sample));
    emf.distribution_fit = &df_mock.get();
    e.improve_fit(ftsr, btsr, nbtsr, edmans, probability, &emf);
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
                   .Using(Close(1.109, TOL),
                          0,
                          Close(1.64 * 2.64 / 3.14159, TOL)))
            .Exactly(1);
    Verify(Method(df_mock, add_sample)
                   .Using(Close(1.109, TOL),
                          1,
                          Close(1.25 * 2.25 / 3.14159, TOL)))
            .Exactly(1);
    VerifyNoOtherInvocations(df_mock);
    emf.distribution_fit = original_dist_fit;  // This avoids breaking cleanup.
}

BOOST_AUTO_TEST_SUITE_END()  // emission_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
