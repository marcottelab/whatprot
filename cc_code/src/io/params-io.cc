/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "params-io.h"

// Standard C++ library headers:
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

// Local project headers:
#include "parameterization/model/channel-model.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using std::ifstream;
using std::ofstream;
using std::setprecision;
using std::string;
using std::vector;
}  // namespace

void write_params(const string& filename,
                  unsigned int num_channels,
                  const vector<SequencingModel>& models,
                  const vector<double>& log_ls) {
    ofstream f(filename);
    f << "p_edman_failure,p_detach,p_initial_block,p_cyclic_block";
    for (unsigned int c = 0; c < num_channels; c++) {
        f << ",ch" << c << ":p_bleach,ch" << c << ":p_dud";
    }
    f << ",log(L)\n";
    f << setprecision(17);
    for (unsigned int i = 0; i < models.size(); i++) {
        f << models[i].p_edman_failure << "," << models[i].p_detach << ","
          << models[i].p_initial_block << "," << models[i].p_cyclic_block;
        for (unsigned int c = 0; c < num_channels; c++) {
            f << "," << models[i].channel_models[c]->p_bleach << ","
              << models[i].channel_models[c]->p_dud;
        }
        f << "," << log_ls[i] << "\n";
    }
    f.close();
}

}  // namespace whatprot
