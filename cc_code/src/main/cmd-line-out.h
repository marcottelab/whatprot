/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_MAIN_CMD_LINE_OUT_H
#define WHATPROT_MAIN_CMD_LINE_OUT_H

namespace whatprot {

void print_bad_inputs();
void print_built_classifier(double time);
void print_finished_basic_setup(double time);
void print_finished_classification(double time);
void print_finished_deduping_dye_tracks(int num, double time);
void print_finished_generating_dye_tracks(int num, double time);
void print_finished_generating_radiometries(int num, double time);
void print_finished_saving_results(double time);
void print_invalid_classifier();
void print_invalid_command();
void print_omp_info();
void print_read_dye_seqs(int num, double time);
void print_read_dye_tracks(int num, double time);
void print_read_radiometries(int num, double time);
void print_total_time(double time);
void print_wrong_number_of_inputs();

}  // namespace whatprot

#endif  // WHATPROT_MAIN_CMD_LINE_OUT_H
