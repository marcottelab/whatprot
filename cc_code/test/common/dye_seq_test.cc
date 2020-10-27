// Author: Matthew Beauregard Smith
#include <boost/test/unit_test.hpp>

#include "common/dye_seq.h"

namespace fluoroseq {

BOOST_AUTO_TEST_SUITE(common_suite);
BOOST_AUTO_TEST_SUITE(dye_seq_suite);

BOOST_AUTO_TEST_CASE(constructor_test) {
    int num_channels = 2;
    string s = "0";
    DyeSeq ds(num_channels, s);
    BOOST_TEST(ds.length == 1);
    BOOST_TEST(ds.num_channels == 2);
    BOOST_TEST(ds.seq != (void *) NULL);
    BOOST_TEST(ds.seq[0] == 0);
}

BOOST_AUTO_TEST_CASE(constructor_empty_string_test) {
    int num_channels = 2;
    string s = "";
    DyeSeq ds(num_channels, s);
    BOOST_TEST(ds.length == 0);
    BOOST_TEST(ds.num_channels == 2);
}

BOOST_AUTO_TEST_CASE(constructor_trailing_dots_test) {
    int num_channels = 2;
    string s = "0..";
    DyeSeq ds(num_channels, s);
    BOOST_TEST(ds.length == 1);
    BOOST_TEST(ds.num_channels == 2);
    BOOST_TEST(ds.seq != (void *) NULL);
    BOOST_TEST(ds.seq[0] == 0);
}

BOOST_AUTO_TEST_CASE(copy_constructor_test) {
    int num_channels = 2;
    string s = "0";
    DyeSeq ds1(num_channels, s);
    DyeSeq ds2(ds1);
    BOOST_TEST(ds1.length == 1);
    BOOST_TEST(ds1.num_channels == 2);
    BOOST_TEST(ds1.seq != (void *) NULL);
    BOOST_TEST(ds1.seq[0] == 0);
    BOOST_TEST(ds2.length == 1);
    BOOST_TEST(ds2.num_channels == 2);
    BOOST_TEST(ds2.seq != (void *) NULL);
    BOOST_TEST(ds2.seq[0] == 0);
}

BOOST_AUTO_TEST_CASE(bracket_op_const_test) {
    int num_channels = 2;
    string s = "0";
    DyeSeq ds(num_channels, s);
    const DyeSeq& cds = ds;
    BOOST_TEST(cds[0] == 0);
}

BOOST_AUTO_TEST_CASE(bracket_op_const_past_end_test) {
    int num_channels = 2;
    string s = "0";
    DyeSeq ds(num_channels, s);
    const DyeSeq& cds = ds;
    BOOST_TEST(cds[1] == -1);
}

BOOST_AUTO_TEST_CASE(more_difficult_dye_seq_string_test) {
    int num_channels = 3;
    string s = "..0.1..2..";
    DyeSeq ds(num_channels, s);
    const DyeSeq& cds = ds;
    BOOST_TEST(cds[0] == -1);
    BOOST_TEST(cds[1] == -1);
    BOOST_TEST(cds[2] == 0);
    BOOST_TEST(cds[3] == -1);
    BOOST_TEST(cds[4] == 1);
    BOOST_TEST(cds[5] == -1);
    BOOST_TEST(cds[6] == -1);
    BOOST_TEST(cds[7] == 2);
    BOOST_TEST(cds[8] == -1);
    BOOST_TEST(cds[9] == -1);
    BOOST_TEST(cds[10] == -1);
}

BOOST_AUTO_TEST_CASE(bracket_op_test) {
    int num_channels = 2;
    string s = "0";
    DyeSeq ds(num_channels, s);
    BOOST_TEST(ds[0] == 0);
    ds[0] = 1;
    BOOST_TEST(ds[0] == 1);
}

BOOST_AUTO_TEST_SUITE_END();  // dye_seq_suite
BOOST_AUTO_TEST_SUITE_END();  // common_suite

}  // namespace fluoroseq
