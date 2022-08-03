/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "scores-io.h"

// Standard C++ library headers:
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <string>
#include <vector>

namespace whatprot {

namespace {
using std::ofstream;
using std::setprecision;
using std::string;
using std::vector;
}  // namespace

void write_scores(const string& filename,
                  const vector<vector<double>>& all_scores) {
    ofstream f(filename);
    for (unsigned int i = 0; i < all_scores.size(); i++) {
        f << setprecision(17) << all_scores[i][0];
        for (unsigned int j = 1; j < all_scores[0].size(); j++) {
            f << "," << setprecision(17) << all_scores[i][j];
        }
        f << "\n";
    }
    f.close();
}

}  // namespace whatprot
