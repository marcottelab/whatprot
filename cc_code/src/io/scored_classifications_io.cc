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

void write_scored_classifications(const string& filename,
                                  int num_radiometries,
                                  const ScoredClassification* results) {
    
    ofstream fpred(filename);
    fpred << "radmat_iz,best_pep_iz,best_pep_score\n";
    for (int i = 0; i < num_radiometries; i++) {
        fpred << i << ",";
        fpred << results[i].id << ",";
        fpred << setprecision(17) << results[i].adjusted_score() << "\n";
    }
    fpred.close();
}

}  // namespace fluoroseq
