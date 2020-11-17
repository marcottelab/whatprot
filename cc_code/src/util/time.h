/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_UTIL_TIME_H
#define FLUOROSEQ_UTIL_TIME_H

namespace fluoroseq {

double wall_time();
double wall_tick();
unsigned int time_based_seed();

}  // namespace fluoroseq

#endif  // FLUOROSEQ_UTIL_TIME_H