// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_MAIN_CMD_LINE_OUT_H
#define FLUOROSEQ_MAIN_CMD_LINE_OUT_H

namespace fluoroseq {

void print_bad_inputs();
void print_built_classifier(double time);
void print_finished_basic_setup(double time);
void print_finished_classification(double time);
void print_finished_saving_results(double time);
void print_invalid_classifier();
void print_mpi_info();
void print_read_dye_seqs(int num, double time);
void print_read_dye_tracks(int num, double time);
void print_read_radiometries(int num, double time);
void print_total_time(double time);
void print_wrong_number_of_inputs();

}  // namespace fluoroseq

#endif  // FLUOROSEQ_MAIN_CMD_LINE_OUT_H
