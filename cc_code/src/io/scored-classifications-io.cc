/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "scored-classifications-io.h"

// Standard C++ library headers:
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <string>
#include <vector>

// Local project headers:
#include "common/scored-classification.h"

namespace whatprot {

namespace {
using std::ofstream;
using std::setprecision;
using std::string;
using std::vector;
}  // namespace

void write_scored_classifications(
        const string& filename,
        int total_num_scored_classifications,
        const vector<ScoredClassification>& scored_classifications) {
    int num_scored_classifications = scored_classifications.size();
    int* ids;
    double* scores;
    convert_raw_from_scored_classifications(
            scored_classifications, &ids, &scores);
    write_scored_classifications_raw(
            filename, total_num_scored_classifications, ids, scores);
}

void convert_raw_from_scored_classifications(
        const vector<ScoredClassification>& scored_classifications,
        int** ids,
        double** scores) {
    *ids = new int[scored_classifications.size()];
    *scores = new double[scored_classifications.size()];
    for (int i = 0; i < scored_classifications.size(); i++) {
        (*ids)[i] = scored_classifications[i].id;
        (*scores)[i] = scored_classifications[i].adjusted_score();
    }
}

void write_scored_classifications_raw(const string& filename,
                                      int num_scored_classifications,
                                      int* ids,
                                      double* scores) {
    ofstream f(filename);
    f << "radmat_iz,best_pep_iz,best_pep_score\n";
    for (int i = 0; i < num_scored_classifications; i++) {
        f << i << ",";
        f << ids[i] << ",";
        f << setprecision(17) << scores[i] << "\n";
    }
    f.close();
    delete[] ids;
    delete[] scores;
}

}  // namespace whatprot
