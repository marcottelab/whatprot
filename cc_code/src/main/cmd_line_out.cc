// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "cmd_line_out.h"

#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

namespace fluoroseq {

namespace {
using std::cout;
}  // namespace

void print_bad_inputs() {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Bad inputs.\n";
}

void print_built_classifier(double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Built classifier (" << time << " seconds).\n";
}

void print_finished_basic_setup(double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Finished basic setup (" << time << " seconds).\n";
}

void print_finished_classification(double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Finished classification (" << time << " seconds).\n";
}

void print_finished_generating_radiometries(int num, double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Finished generating " << num << " radiometries (" << time
         << "seconds).\n";
}

void print_finished_saving_results(double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Finished saving results (" << time << " seconds).\n";
}

void print_invalid_classifier() {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Invalid classifier. First argument must be 'hmm', 'ann', or "
         << "'hybrid'.\n";
}

void print_mpi_info() {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    cout << "Using MPI with " << mpi_size << " processes.\n";
#else  // USE_MPI
    cout << "Not using MPI.\n";
#endif  // USE_MPI
}

void print_read_dye_seqs(int num, double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Read " << num << " dye seqs (" << time << " seconds).\n";
}

void print_read_dye_tracks(int num, double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Read " << num << " dye tracks (" << time << " seconds).\n";
}

void print_read_radiometries(int num, double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Read " << num << " radiometries (" << time << " seconds).\n";
}

void print_total_time(double time) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Total run time: " << time << " seconds.\n";
}

void print_wrong_number_of_inputs() {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    cout << "Wrong number of inputs.\n";
}

}  // namespace fluoroseq
