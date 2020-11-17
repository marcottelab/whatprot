/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// This file creates the main function for the boost testing framework. The
// boost testing framework can only be included in this manner one time - all
// other inclusions are of <boost/test/unit_test.hpp> which is NOT the same.

// Declare module first:
#define BOOST_TEST_MODULE fluoroseq

// Then include:
#include <boost/test/included/unit_test.hpp>
