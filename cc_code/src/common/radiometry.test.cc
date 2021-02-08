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
#include "radiometry.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(common_suite)
BOOST_AUTO_TEST_SUITE(radiometry_suite)

BOOST_AUTO_TEST_CASE(constructor_test) {
    int num_timesteps = 3;
    int num_channels = 2;
    Radiometry r(num_timesteps, num_channels);
    BOOST_TEST(r.num_timesteps = num_timesteps);
    BOOST_TEST(r.num_channels = num_channels);
    BOOST_TEST(r.intensities != (void*)NULL);
}

BOOST_AUTO_TEST_CASE(copy_constructor_test) {
    int num_timesteps = 3;
    int num_channels = 2;
    Radiometry r1(num_timesteps, num_channels);
    r1.intensities[0] = 60;
    r1.intensities[1] = 61;
    r1.intensities[2] = 62;
    r1.intensities[3] = 63;
    r1.intensities[4] = 64;
    r1.intensities[5] = 65;
    Radiometry r2(r1);
    BOOST_TEST(r1.num_timesteps = num_timesteps);
    BOOST_TEST(r1.num_channels = num_channels);
    BOOST_TEST(r1.intensities != (void*)NULL);
    BOOST_TEST(r1.intensities[0] == 60);
    BOOST_TEST(r1.intensities[1] == 61);
    BOOST_TEST(r1.intensities[2] == 62);
    BOOST_TEST(r1.intensities[3] == 63);
    BOOST_TEST(r1.intensities[4] == 64);
    BOOST_TEST(r1.intensities[5] == 65);
    BOOST_TEST(r2.num_timesteps = num_timesteps);
    BOOST_TEST(r2.num_channels = num_channels);
    BOOST_TEST(r2.intensities != (void*)NULL);
    BOOST_TEST(r2.intensities[0] == 60);
    BOOST_TEST(r2.intensities[1] == 61);
    BOOST_TEST(r2.intensities[2] == 62);
    BOOST_TEST(r2.intensities[3] == 63);
    BOOST_TEST(r2.intensities[4] == 64);
    BOOST_TEST(r2.intensities[5] == 65);
}

BOOST_AUTO_TEST_CASE(paren_op_test) {
    int num_timesteps = 3;
    int num_channels = 2;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 600;
    r(0, 1) = 601;
    r(1, 0) = 610;
    r(1, 1) = 611;
    r(2, 0) = 620;
    r(2, 1) = 621;
    BOOST_TEST(r(0, 0) == 600);
    BOOST_TEST(r(0, 1) == 601);
    BOOST_TEST(r(1, 0) == 610);
    BOOST_TEST(r(1, 1) == 611);
    BOOST_TEST(r(2, 0) == 620);
    BOOST_TEST(r(2, 1) == 621);
}

BOOST_AUTO_TEST_CASE(paren_op_const_test) {
    int num_timesteps = 3;
    int num_channels = 2;
    Radiometry r(num_timesteps, num_channels);
    r(0, 0) = 600;
    r(0, 1) = 601;
    r(1, 0) = 610;
    r(1, 1) = 611;
    r(2, 0) = 620;
    r(2, 1) = 621;
    const Radiometry& cr = r;
    BOOST_TEST(cr(0, 0) == 600);
    BOOST_TEST(cr(0, 1) == 601);
    BOOST_TEST(cr(1, 0) == 610);
    BOOST_TEST(cr(1, 1) == 611);
    BOOST_TEST(cr(2, 0) == 620);
    BOOST_TEST(cr(2, 1) == 621);
}

BOOST_AUTO_TEST_SUITE_END()  // radiometry_suite
BOOST_AUTO_TEST_SUITE_END()  // common_suite

}  // namespace whatprot
