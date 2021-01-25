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

// Local project headers:
#include "common/radiometry.h"
#include "tensor/tensor.h"

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
using std::function;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
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

BOOST_AUTO_TEST_SUITE_END()  // emission_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace fluoroseq
