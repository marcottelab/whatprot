// Author: Matthew Beauregard Smith (UT Austin)
#include "scored_classifications_io.h"

#include <fstream>
#include <iomanip>  // for std::setprecision
#include <string>

#include "common/scored_classification.h"

namespace fluoroseq {

namespace {
using std::ofstream;
using std::setprecision;
using std::string;
}  // namespace

void write_scored_classifications(
        const string& filename,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications) {
    ofstream f(filename);
    f << "radmat_iz,best_pep_iz,best_pep_score\n";
    for (int i = 0; i < num_scored_classifications; i++) {
        f << i << ",";
        f << scored_classifications[i].id << ",";
        f << setprecision(17)
          << scored_classifications[i].adjusted_score()
          << "\n";
    }
    f.close();
}

}  // namespace fluoroseq
