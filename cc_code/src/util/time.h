/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_UTIL_TIME_H
#define WHATPROT_UTIL_TIME_H

namespace whatprot {

double wall_time();
double wall_tick();
unsigned int time_based_seed();

}  // namespace whatprot

#endif  // WHATPROT_UTIL_TIME_H