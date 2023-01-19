/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// The tests in this file may need to be run with valgrind to detect errors.

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// File under test:
#include "schroedinger-pointers.h"

namespace whatprot {

BOOST_AUTO_TEST_SUITE(util_suite)
BOOST_AUTO_TEST_SUITE(delete_suite)

BOOST_AUTO_TEST_CASE(dereference_if_pointer_not_pointer_test) {
    int x = 42;
    int y = dereference_if_pointer(x);
    BOOST_TEST(y == 42);
}

BOOST_AUTO_TEST_CASE(dereference_if_pointer_is_pointer_test) {
    int x = 42;
    int* y = &x;
    int z = dereference_if_pointer(y);
    BOOST_TEST(z == 42);
}

BOOST_AUTO_TEST_CASE(delete_if_pointer_not_pointer_test) {
    int x = 42;
    delete_if_pointer(x);
}

BOOST_AUTO_TEST_CASE(delete_if_pointer_is_pointer_test) {
    int* x = new int;
    *x = 42;
    delete_if_pointer(x);
}

BOOST_AUTO_TEST_CASE(delete_array_not_pointers_test) {
    unsigned int size = 3;
    int* arr = new int[size];
    delete_array(size, arr);
}

BOOST_AUTO_TEST_CASE(delete_array_is_pointers_test) {
    unsigned int size = 3;
    int** arr = new int*[size];
    for (unsigned int i = 0; i < size; i++) {
        arr[i] = new int;
        *arr[i] = i;
    }
    delete_array(size, arr);
}

BOOST_AUTO_TEST_SUITE_END()  // delete_suite
BOOST_AUTO_TEST_SUITE_END()  // util_suite

}  // namespace whatprot
