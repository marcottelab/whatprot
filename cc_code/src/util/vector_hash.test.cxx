// Author: Matthew Beauregard Smith
//
// Yes these tests can fail for valid reasons, but it is highly unlikely. These
// tests still serve as a good sanity check that things are working correctly.
#include <boost/test/unit_test.hpp>

#include "util/vector_hash.h"

#include <functional>
#include <vector>

namespace fluoroseq {

namespace {
using std::vector;
using std::hash;
}

BOOST_AUTO_TEST_SUITE(util_suite);
BOOST_AUTO_TEST_SUITE(vector_hash_suite);

BOOST_AUTO_TEST_CASE(same_num_same_hash_test) {
    vector<int> v1;
    v1.push_back(42);
    vector<int> v2;
    v2.push_back(42);
    hash<vector<int>> hasher;
    BOOST_TEST(hasher(v1) == hasher(v2));
}

BOOST_AUTO_TEST_CASE(same_larger_vectors_same_hash_test) {
    vector<int> v1;
    v1.push_back(42);
    v1.push_back(13);
    v1.push_back(7);
    vector<int> v2;
    v2.push_back(42);
    v2.push_back(13);
    v2.push_back(7);
    hash<vector<int>> hasher;
    BOOST_TEST(hasher(v1) == hasher(v2));
}

BOOST_AUTO_TEST_CASE(different_num_different_hash_test) {
    vector<int> v1;
    v1.push_back(1);
    vector<int> v2;
    v2.push_back(2);
    hash<vector<int>> hasher;
    BOOST_TEST(hasher(v1) != hasher(v2));
}

BOOST_AUTO_TEST_CASE(different_size_different_hash_test) {
    vector<int> v1;
    v1.push_back(42);
    vector<int> v2;
    v2.push_back(42);
    v2.push_back(42);
    hash<vector<int>> hasher;
    BOOST_TEST(hasher(v1) != hasher(v2));
}

BOOST_AUTO_TEST_CASE(different_element_order_different_hash_test) {
    vector<int> v1;
    v1.push_back(42);
    v1.push_back(13);
    v1.push_back(7);
    vector<int> v2;
    v2.push_back(13);
    v2.push_back(42);
    v2.push_back(7);
    hash<vector<int>> hasher;
    BOOST_TEST(hasher(v1) != hasher(v2));
}

BOOST_AUTO_TEST_SUITE_END();  // vector_hash_suite
BOOST_AUTO_TEST_SUITE_END();  // util_suite

}  // namespace fluoroseq
