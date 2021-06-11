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

namespace whatprot {

namespace {
using std::ifstream;
using std::ofstream;
using std::setprecision;
using std::string;
using std::vector;
}  // namespace

void read_radiometries(const string& filename,
                       unsigned int* num_timesteps,
                       unsigned int* num_channels,
                       unsigned int* total_num_radiometries,
                       vector<Radiometry>* radiometries) {
    unsigned int num_radiometries;
    double* intensities;
    read_radiometries_raw(filename,
                          num_timesteps,
                          num_channels,
                          total_num_radiometries,
                          &num_radiometries,
                          &intensities);
    num_radiometries = *total_num_radiometries;
    convert_radiometries_from_raw(*num_timesteps,
                                  *num_channels,
                                  num_radiometries,
                                  intensities,
                                  radiometries);
    delete[] intensities;
}

void read_radiometries_raw(const string& filename,
                           unsigned int* num_timesteps,
                           unsigned int* num_channels,
                           unsigned int* total_num_radiometries,
                           unsigned int* num_radiometries,
                           double** intensities) {
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> *total_num_radiometries;
    *num_radiometries = *total_num_radiometries;
    *intensities = new double[(*total_num_radiometries) * (*num_timesteps)
                              * (*num_channels)];
    for (unsigned int i = 0;
         i < (*total_num_radiometries) * (*num_timesteps) * (*num_channels);
         i++) {
        f >> (*intensities)[i];
    }
    f.close();
}

void convert_radiometries_from_raw(unsigned int num_timesteps,
                                   unsigned int num_channels,
                                   unsigned int num_radiometries,
                                   double* intensities,
                                   vector<Radiometry>* radiometries) {
    radiometries->reserve(num_radiometries);
    for (unsigned int i = 0; i < num_radiometries; i++) {
        radiometries->push_back(Radiometry(num_timesteps, num_channels));
        for (unsigned int j = 0; j < num_timesteps * num_channels; j++) {
            radiometries->back().intensities[j] =
                    intensities[i * (num_timesteps * num_channels) + j];
        }
    }
}

void write_radiometries(
        const string& filename,
        unsigned int num_timesteps,
        unsigned int num_channels,
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries) {
    double* intensities;
    convert_raw_from_radiometries(
            radiometries,
            num_timesteps * num_channels,  // radiometry size
            &intensities);
    int total_num_radiometries = radiometries.size();
    write_radiometries_raw(filename,
                           num_timesteps,
                           num_channels,
                           total_num_radiometries,
                           intensities);
}

void convert_raw_from_radiometries(
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries,
        unsigned int radiometry_size,
        double** intensities) {
    *intensities = new double[radiometries.size() * radiometry_size];
    for (unsigned int i = 0; i < radiometries.size(); i++) {
        for (unsigned int j = 0; j < radiometry_size; j++) {
            (*intensities)[i * radiometry_size + j] =
                    radiometries[i].value.intensities[j];
        }
    }
}

void write_radiometries_raw(const std::string& filename,
                            unsigned int num_timesteps,
                            unsigned int num_channels,
                            unsigned int num_radiometries,
                            double* intensities) {
    ofstream f(filename);
    f << num_timesteps << "\n";
    f << num_channels << "\n";
    f << num_radiometries << "\n";
    for (unsigned int i = 0; i < num_radiometries; i++) {
        for (unsigned int j = 0; j < num_timesteps * num_channels; j++) {
            if (j != 0) {
                f << "\t";
            }
            f << setprecision(17)
              << intensities[i * num_timesteps * num_channels + j];
        }
        f << "\n";
    }
    f.close();
    delete[] intensities;
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
