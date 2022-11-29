/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "radiometries-io.h"

// Standard C++ library headers:
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <string>
#include <vector>

// Local project headers:
#include "common/radiometry.h"
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

void read_radiometries(const string& filename,
                       const SequencingModel& seq_model,
                       unsigned int* num_timesteps,
                       unsigned int* num_channels,
                       unsigned int* num_radiometries,
                       vector<Radiometry>* radiometries) {
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> *num_radiometries;
    radiometries->reserve(*num_radiometries);
    for (unsigned int i = 0; i < (*num_radiometries); i++) {
        radiometries->push_back(Radiometry(*num_timesteps, *num_channels));
        for (unsigned int t = 0; t < (*num_timesteps); t++) {
            for (unsigned int c = 0; c < (*num_channels); c++) {
                double intensity;
                f >> intensity;
                radiometries->back()(t, c) =
                        intensity / seq_model.channel_models[c]->mu;
            }
        }
    }
    f.close();
}

void write_radiometries(
        const string& filename,
        const SequencingModel& seq_model,
        unsigned int num_timesteps,
        unsigned int num_channels,
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries) {
    ofstream f(filename);
    f << num_timesteps << "\n";
    f << num_channels << "\n";
    f << radiometries.size() << "\n";
    for (unsigned int i = 0; i < radiometries.size(); i++) {
        for (unsigned int t = 0; t < num_timesteps; t++) {
            for (unsigned int c = 0; c < num_channels; c++) {
                if (t != 0 || c != 0) {
                    f << "\t";
                }
                f << setprecision(17)
                  << seq_model.channel_models[c]->mu
                                * radiometries[i].value(t, c);
            }
        }
        f << "\n";
    }
    f.close();
}

void write_ys(
        const string& filename,
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries) {
    int* ys;
    get_raw_ys(radiometries, &ys);
    int total_num_radiometries = radiometries.size();
    write_ys_raw(filename,
                 total_num_radiometries,  // num radiometries
                 ys);
}

void get_raw_ys(
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries,
        int** ys) {
    *ys = new int[radiometries.size()];
    for (unsigned int i = 0; i < radiometries.size(); i++) {
        (*ys)[i] = radiometries[i].source.source;
    }
}

void write_ys_raw(const string& filename,
                  unsigned int num_radiometries,
                  int* ys) {
    ofstream f(filename);
    f << num_radiometries << "\n";
    for (unsigned int i = 0; i < num_radiometries; i++) {
        f << ys[i] << "\n";
    }
    f.close();
    delete[] ys;
}

}  // namespace whatprot
