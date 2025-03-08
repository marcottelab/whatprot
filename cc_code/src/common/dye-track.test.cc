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
#include "dye-track.h"

// Standard C++ library headers:
#include <utility>

// Local project headers:
#include "common/dye-seq.h"

namespace whatprot {

namespace {
using std::hash;
using std::move;
}  // namespace

BOOST_AUTO_TEST_SUITE(common_suite)
BOOST_AUTO_TEST_SUITE(dye_track_suite)

BOOST_AUTO_TEST_CASE(constructor_from_dye_seq_trivial_test) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    DyeSeq ds(num_channels, "0");
    DyeTrack dt(num_timesteps, num_channels, ds);
    BOOST_TEST(dt.num_timesteps == num_timesteps);
    BOOST_TEST(dt.num_channels == num_channels);
    BOOST_TEST(dt.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt(0, 0) == 1);
}

BOOST_AUTO_TEST_CASE(constructor_from_dye_seq_bigger_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01..1");
    DyeTrack dt(num_timesteps, num_channels, ds);
    BOOST_TEST(dt.num_timesteps == num_timesteps);
    BOOST_TEST(dt.num_channels == num_channels);
    BOOST_TEST(dt.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt(0, 0) == 1);
    BOOST_TEST(dt(0, 1) == 2);
    BOOST_TEST(dt(1, 0) == 0);
    BOOST_TEST(dt(1, 1) == 2);
    BOOST_TEST(dt(2, 0) == 0);
    BOOST_TEST(dt(2, 1) == 1);
}

BOOST_AUTO_TEST_CASE(constructor_from_dye_seq_many_timesteps_few_aas_test) {
    unsigned int num_timesteps = 5;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01");
    DyeTrack dt(num_timesteps, num_channels, ds);
    BOOST_TEST(dt.num_timesteps == num_timesteps);
    BOOST_TEST(dt.num_channels == num_channels);
    BOOST_TEST(dt.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt(0, 0) == 1);
    BOOST_TEST(dt(0, 1) == 1);
    BOOST_TEST(dt(1, 0) == 0);
    BOOST_TEST(dt(1, 1) == 1);
    BOOST_TEST(dt(2, 0) == 0);
    BOOST_TEST(dt(2, 1) == 0);
    BOOST_TEST(dt(3, 0) == 0);
    BOOST_TEST(dt(3, 1) == 0);
    BOOST_TEST(dt(4, 0) == 0);
    BOOST_TEST(dt(4, 1) == 0);
}

BOOST_AUTO_TEST_CASE(constructor_from_short_array_trivial_test) {
    unsigned int num_timesteps = 1;
    unsigned int num_channels = 1;
    short* counts = new short[num_timesteps * num_channels];
    counts[0] = 1;
    DyeTrack dt(num_timesteps, num_channels, counts);
    delete[] counts;
    BOOST_TEST(dt.num_timesteps == num_timesteps);
    BOOST_TEST(dt.num_channels == num_channels);
    BOOST_TEST(dt.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt(0, 0) == 1);
}

BOOST_AUTO_TEST_CASE(constructor_from_short_array_bigger_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    short* counts = new short[num_timesteps * num_channels];
    counts[0] = 1;
    counts[1] = 2;
    counts[2] = 0;
    counts[3] = 2;
    counts[4] = 0;
    counts[5] = 1;
    DyeTrack dt(num_timesteps, num_channels, counts);
    delete[] counts;
    BOOST_TEST(dt.num_timesteps == num_timesteps);
    BOOST_TEST(dt.num_channels == num_channels);
    BOOST_TEST(dt.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt(0, 0) == 1);
    BOOST_TEST(dt(0, 1) == 2);
    BOOST_TEST(dt(1, 0) == 0);
    BOOST_TEST(dt(1, 1) == 2);
    BOOST_TEST(dt(2, 0) == 0);
    BOOST_TEST(dt(2, 1) == 1);
}

BOOST_AUTO_TEST_CASE(constructor_no_data_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeTrack dt(num_timesteps, num_channels);
    BOOST_TEST(dt.num_timesteps == num_timesteps);
    BOOST_TEST(dt.num_channels == num_channels);
    BOOST_TEST(dt.counts.size() == num_timesteps * num_channels);
}

BOOST_AUTO_TEST_CASE(copy_constructor_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01..1");
    DyeTrack dt1(num_timesteps, num_channels, ds);
    DyeTrack dt2(dt1);
    BOOST_TEST(dt1.num_timesteps == num_timesteps);
    BOOST_TEST(dt1.num_channels == num_channels);
    BOOST_TEST(dt1.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt1(0, 0) == 1);
    BOOST_TEST(dt1(0, 1) == 2);
    BOOST_TEST(dt1(1, 0) == 0);
    BOOST_TEST(dt1(1, 1) == 2);
    BOOST_TEST(dt1(2, 0) == 0);
    BOOST_TEST(dt1(2, 1) == 1);
    BOOST_TEST(dt2.num_timesteps == num_timesteps);
    BOOST_TEST(dt2.num_channels == num_channels);
    BOOST_TEST(dt2.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt2(0, 0) == 1);
    BOOST_TEST(dt2(0, 1) == 2);
    BOOST_TEST(dt2(1, 0) == 0);
    BOOST_TEST(dt2(1, 1) == 2);
    BOOST_TEST(dt2(2, 0) == 0);
    BOOST_TEST(dt2(2, 1) == 1);
}

BOOST_AUTO_TEST_CASE(move_constructor_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01..1");
    DyeTrack dt1(num_timesteps, num_channels, ds);
    DyeTrack dt2(move(dt1));
    BOOST_TEST(dt2.num_timesteps == num_timesteps);
    BOOST_TEST(dt2.num_channels == num_channels);
    BOOST_TEST(dt2.counts.size() == num_timesteps * num_channels);
    BOOST_TEST(dt2(0, 0) == 1);
    BOOST_TEST(dt2(0, 1) == 2);
    BOOST_TEST(dt2(1, 0) == 0);
    BOOST_TEST(dt2(1, 1) == 2);
    BOOST_TEST(dt2(2, 0) == 0);
    BOOST_TEST(dt2(2, 1) == 1);
}

BOOST_AUTO_TEST_CASE(equal_op_true_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01..1");
    DyeTrack dt1(num_timesteps, num_channels, ds);
    DyeTrack dt2(num_timesteps, num_channels, ds);
    BOOST_TEST((dt1 == dt2));
}

BOOST_AUTO_TEST_CASE(equal_op_false_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeSeq ds1(num_channels, "01..1");
    DyeSeq ds2(num_channels, "0.0..1");
    DyeTrack dt1(num_timesteps, num_channels, ds1);
    DyeTrack dt2(num_timesteps, num_channels, ds2);
    BOOST_TEST(!(dt1 == dt2));
}

BOOST_AUTO_TEST_CASE(paren_op_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeTrack dt(num_timesteps, num_channels);
    dt(1, 1) = 4;
    BOOST_TEST(dt(1, 1) == 4);
}

BOOST_AUTO_TEST_CASE(paren_op_const_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeTrack dt(num_timesteps, num_channels);
    dt(1, 1) = 4;
    const DyeTrack& cdt = dt;
    BOOST_TEST(cdt(1, 1) == 4);
}

BOOST_AUTO_TEST_CASE(bracket_op_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeTrack dt(num_timesteps, num_channels);
    dt[1][1] = 4;
    BOOST_TEST(dt[1][1] == 4);
}

BOOST_AUTO_TEST_CASE(bracket_op_const_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeTrack dt(num_timesteps, num_channels);
    dt[1][1] = 4;
    const DyeTrack& cdt = dt;
    BOOST_TEST(cdt[1][1] == 4);
}

// This test may fail by bad luck, but this is very unlikely.
BOOST_AUTO_TEST_CASE(hash_same_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeSeq ds(num_channels, "01..1");
    DyeTrack dt1(num_timesteps, num_channels, ds);
    DyeTrack dt2(num_timesteps, num_channels, ds);
    hash<DyeTrack> hasher;
    BOOST_TEST(hasher(dt1) == hasher(dt2));
}

// This test may fail by bad luck, but this is very unlikely.
BOOST_AUTO_TEST_CASE(hash_different_test) {
    unsigned int num_timesteps = 3;
    unsigned int num_channels = 2;
    DyeSeq ds1(num_channels, "01..1");
    DyeSeq ds2(num_channels, "0.0..1");
    DyeTrack dt1(num_timesteps, num_channels, ds1);
    DyeTrack dt2(num_timesteps, num_channels, ds2);
    hash<DyeTrack> hasher;
    BOOST_TEST(hasher(dt1) != hasher(dt2));
}

BOOST_AUTO_TEST_SUITE_END()  // dye_track_suite
BOOST_AUTO_TEST_SUITE_END()  // common_suite

}  // namespace whatprot
