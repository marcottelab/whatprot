/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "cmd-line-out.h"

// Standard C++ library headers:
#include <iostream>

// OpenMP
#include <omp.h>

namespace whatprot {

namespace {
using std::cout;
}  // namespace

void print_bad_inputs() {
    cout << "Bad inputs.\n";
}

void print_built_classifier(double time) {
    cout << "Built classifier (" << time << " seconds).\n";
}

void print_finished_basic_setup(double time) {
    cout << "Finished basic setup (" << time << " seconds).\n";
}

void print_finished_classification(double time) {
    cout << "Finished classification (" << time << " seconds).\n";
}

void print_finished_deduping_dye_tracks(int num, double time) {
    cout << "Finished deduping " << num << " unique dye tracks (" << time
         << "seconds).\n";
}

void print_finished_generating_dye_tracks(int num, double time) {
    int total_num = num;
    cout << "Finished generating " << total_num << " dye tracks (" << time
         << "seconds).\n";
}

void print_finished_generating_radiometries(int num, double time) {
    int total_num = num;
    cout << "Finished generating " << total_num << " radiometries (" << time
         << "seconds).\n";
}

void print_finished_saving_results(double time) {
    cout << "Finished saving results (" << time << " seconds).\n";
}

void print_invalid_classifier() {
    cout << "Invalid classifier. Second argument must be 'hmm', 'nn', or "
         << "'hybrid'.\n";
}

void print_invalid_command() {
    cout << "Invalid command. First argument must be 'classify' or "
         << "'simulate'.\n";
}

void print_omp_info() {
    cout << "Using OpenMP with " << omp_get_max_threads() << " threads.\n";
}

void print_read_dye_seqs(int num, double time) {
    cout << "Read " << num << " dye seqs (" << time << " seconds).\n";
}

void print_read_dye_tracks(int num, double time) {
    cout << "Read " << num << " dye tracks (" << time << " seconds).\n";
}

void print_read_radiometries(int num, double time) {
    cout << "Read " << num << " radiometries (" << time << " seconds).\n";
}

void print_total_time(double time) {
    cout << "Total run time: " << time << " seconds.\n";
}

void print_wrong_number_of_inputs() {
    cout << "Wrong number of inputs.\n";
}

}  // namespace whatprot
