// Author: Matthew Beauregard Smith (UT Austin)
#include "scored_classifications_mpi_io.h"

#include <fstream>
#include <iomanip>  // for std::setprecision
#include <string>

#include <mpi.h>

#include "common/scored_classification.h"

namespace fluoroseq {

namespace {
using std::ofstream;
using std::setprecision;
using std::string;
}  // namespace

void mpi_write_scored_classifications(
        const std::string& filename,
        int total_num_scored_classifications,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == 0) {
        mpi_write_scored_classifications_master(
                filename,
                total_num_scored_classifications,
                num_scored_classifications,
                scored_classifications);
    } else {
        mpi_write_scored_classifications_slave(
                filename,
                total_num_scored_classifications,
                num_scored_classifications,
                scored_classifications);
    }
}

void mpi_write_scored_classifications_master(
        const std::string& filename,
        int total_num_scored_classifications,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications) {
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    int* all_ids = new int[total_num_scored_classifications];
    double* all_scores = new double[total_num_scored_classifications];
    int* ids = new int[num_scored_classifications];
    double* scores = new double[num_scored_classifications];
    for (int i = 0; i < num_scored_classifications; i++) {
        ids[i] = scored_classifications[i].id;
        scores[i] = scored_classifications[i].adjusted_score();
    }
    int* result_recvcounts = new int[mpi_size];
    int* result_displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long) total_num_scored_classifications * (long) i
                    / (long) mpi_size;
        int end = (long) total_num_scored_classifications * (long) (i + 1)
                  / (long) mpi_size;
        result_recvcounts[i] = end - begin;
        result_displs[i] = begin;
    }
    MPI_Gatherv(ids,
                num_scored_classifications,
                MPI_INT,
                all_ids,
                result_recvcounts,
                result_displs,
                MPI_INT,
                0,  // root
                MPI_COMM_WORLD);
    MPI_Gatherv(scores,
                num_scored_classifications,
                MPI_DOUBLE,
                all_scores,
                result_recvcounts,
                result_displs,
                MPI_DOUBLE,
                0,  // root
                MPI_COMM_WORLD);
    ofstream f(filename);
    f << "radmat_iz,best_pep_iz,best_pep_score\n";
    for (int i = 0; i < total_num_scored_classifications; i++) {
        f << i << ",";
        f << all_ids[i] << ",";
        f << setprecision(17)
          << all_scores[i]
          << "\n";
    }
    f.close();
    delete[] all_ids;
    delete[] all_scores;
    delete[] ids;
    delete[] scores;
    delete[] result_recvcounts;
    delete[] result_displs;
}

void mpi_write_scored_classifications_slave(
        const std::string& filename,
        int total_num_scored_classifications,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications) {
    int* ids = new int[num_scored_classifications];
    double* scores = new double[num_scored_classifications];
    for (int i = 0; i < num_scored_classifications; i++) {
        ids[i] = scored_classifications[i].id;
        scores[i] = scored_classifications[i].adjusted_score();
    }
    MPI_Gatherv(ids,
                num_scored_classifications,
                MPI_INT,
                NULL,  // recvbuf
                NULL,  // recvcounts
                NULL,  // displs
                MPI_INT,
                0,  // root
                MPI_COMM_WORLD);
    MPI_Gatherv(scores,
                num_scored_classifications,
                MPI_DOUBLE,
                NULL,  // recvbuf
                NULL,  // recvcounts
                NULL,  // displs
                MPI_DOUBLE,
                0,  // root
                MPI_COMM_WORLD);
    delete[] ids;
    delete[] scores;
}

}  // namespace fluoroseq
