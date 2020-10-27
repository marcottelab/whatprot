// Author: Matthew Beauregard Smith
//
// Don't put too much faith in these tests. The "sourced_data.h" file contains
// a lot of bad code and needs to be ripped out and replaced. To the best of my
// knowledge, none of the broken things are being used by other pieces of code.
// I plan to test that these other pieces of code work as expected, which should
// cover the most important things in "sourced_data.h". Nevertheless, it should
// be a priority to get rid of this file.
#include <boost/test/unit_test.hpp>

#include "common/sourced_data.h"

namespace fluoroseq {

BOOST_AUTO_TEST_SUITE(common_suite);
BOOST_AUTO_TEST_SUITE(sourced_data_suite);

BOOST_AUTO_TEST_CASE(sourced_data_value_not_ptr_source_not_ptr_test) {
    double value = 3.14;
    int source = 42;
    SourcedData<double, int> sd(value, source);
    BOOST_TEST(sd.value == 3.14);
    BOOST_TEST(sd.source == 42);
}

BOOST_AUTO_TEST_CASE(sourced_data_value_is_ptr_source_not_ptr_test) {
    double* value = new double;
    *value = 3.14;
    int source = 42;
    SourcedData<double*, int> sd(value, source);
    BOOST_TEST(*sd.value == 3.14);
    BOOST_TEST(sd.source == 42);
}

BOOST_AUTO_TEST_CASE(sourced_data_value_not_ptr_source_is_ptr_test) {
    double value = 3.14;
    int* source = new int;
    *source = 42;
    SourcedData<double, int*> sd(value, source);
    BOOST_TEST(sd.value = 3.14);
    BOOST_TEST(*sd.source = 42);
}

BOOST_AUTO_TEST_CASE(sourced_data_value_is_ptr_source_is_ptr_test) {
    double* value = new double;
    *value = 3.14;
    int* source = new int;
    *source = 42;
    SourcedData<double*, int*> sd(value, source);
    BOOST_TEST(*sd.value == 3.14);
    BOOST_TEST(*sd.source == 42);
}

BOOST_AUTO_TEST_CASE(source_set_sources_not_ptrs_test) {
    int num_sources = 3;
    int* sources = new int[num_sources];
    sources[0] = 30;
    sources[1] = 31;
    sources[2] = 32;
    SourceSet<int> ss(num_sources, sources);
    BOOST_TEST(ss.num_sources == num_sources);
    BOOST_TEST(ss.sources[0] == 30);
    BOOST_TEST(ss.sources[1] == 31);
    BOOST_TEST(ss.sources[2] == 32);
}

BOOST_AUTO_TEST_CASE(source_set_sources_are_ptrs_test) {
    int num_sources = 3;
    int** sources = new int*[num_sources];
    sources[0] = new int;
    *sources[0] = 30;
    sources[1] = new int;
    *sources[1] = 31;
    sources[2] = new int;
    *sources[2] = 32;
    SourceSet<int*> ss(num_sources, sources);
    BOOST_TEST(ss.num_sources == num_sources);
    BOOST_TEST(*ss.sources[0] == 30);
    BOOST_TEST(*ss.sources[1] == 31);
    BOOST_TEST(*ss.sources[2] == 32);
}

BOOST_AUTO_TEST_CASE(source_count_source_not_ptr_test) {
    int source = 42;
    int count = 88;
    SourceCount<int> sc(source, count);
    BOOST_TEST(sc.source == 42);
    BOOST_TEST(sc.count == 88);
}

BOOST_AUTO_TEST_CASE(source_count_source_is_ptr_test) {
    int* source = new int;
    *source = 42;
    int count = 88;
    SourceCount<int*> sc(source, count);
    BOOST_TEST(*sc.source == 42);
    BOOST_TEST(sc.count == 88);
}

BOOST_AUTO_TEST_SUITE_END();  // sourced_data_suite
BOOST_AUTO_TEST_SUITE_END();  // common_suite

}  // namespace fluoroseq
