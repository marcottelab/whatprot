// Author: Matthew Beauregard Smith (UT Austin)
#include "time.h"

#include <ctime>

namespace fluoroseq {

double wall_time() {
    return (double) clock() / (double) CLOCKS_PER_SEC;
}

}  // namespace fluoroseq
