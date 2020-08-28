// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_UTIL_TIME_H
#define FLUOROSEQ_UTIL_TIME_H

namespace fluoroseq {

double wall_time();
double wall_tick();
unsigned int time_based_seed();

}  // namespace fluoroseq

#endif  // FLUOROSEQ_UTIL_TIME_H