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
#include <string>
#include <vector>

// Local project headers:
#include "parameterization/model/channel-model.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
}  // namespace

void write_params(const string& filename,
                  unsigned int num_channels,
                  const vector<SequencingModel>& models) {
    ofstream f(filename);
    f << "p_edman_failure,p_detach,p_initial_block,p_cyclic_block";
    for (unsigned int c = 0; c < num_channels; c++) {
        f << ",ch" << c << ":p_bleach,ch" << c << ":p_dud";
    }
    f << "\n";
    for (const SequencingModel& model : models) {
        f << model.p_edman_failure << "," << model.p_detach << ","
          << model.p_initial_block << "," << model.p_cyclic_block;
        for (unsigned int c = 0; c < num_channels; c++) {
            f << "," << model.channel_models[c]->p_bleach << ","
              << model.channel_models[c]->p_dud;
        }
        f << "\n";
    }
}

}  // namespace whatprot
