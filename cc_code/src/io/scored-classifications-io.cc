/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// For MPI version, define compiler macro USE_MPI when building.

// Defining symbols from header:
#include "scored-classifications-io.h"

// Standard C++ library headers:
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <string>
#include <vector>

// MPI header:
#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

// Local project headers:
#include "common/scored-classification.h"
#ifdef USE_MPI
#include "io/mpi-counts-displs.h"
#endif  // USE_MPI

namespace fluoroseq {

namespace {
using std::ofstream;
using std::setprecision;
using std::string;
using std::vector;
}  // namespace

void write_scored_classifications(
        const string& filename,
        int total_num_scored_classifications,
        const vector<ScoredClassification>& scored_classifications) {
    int num_scored_classifications = scored_classifications.size();
    int* ids;
    double* scores;
    convert_raw_from_scored_classifications(
            scored_classifications, &ids, &scores);
#ifdef USE_MPI
    gather_scored_classifications(total_num_scored_classifications,
                                  num_scored_classifications,
                                  &ids,
                                  &scores);
#endif  // USE_MPI
    write_scored_classifications_raw(
            filename, total_num_scored_classifications, ids, scores);
}

void convert_raw_from_scored_classifications(
        const vector<ScoredClassification>& scored_classifications,
        int** ids,
        double** scores) {
    *ids = new int[scored_classifications.size()];
    *scores = new double[scored_classifications.size()];
    for (int i = 0; i < scored_classifications.size(); i++) {
        (*ids)[i] = scored_classifications[i].id;
        (*scores)[i] = scored_classifications[i].adjusted_score();
    }
}

#ifdef USE_MPI
void gather_scored_classifications(int total_num_scored_classifications,
                                   int num_scored_classifications,
                                   int** ids,
                                   double** scores) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int* mpi_counts;
    int* mpi_displs;
    mpi_counts_displs(total_num_scored_classifications,
                      1,  // block size
                      &mpi_counts,
                      &mpi_displs);

    int* ids_recv_buf;
    double* scores_recv_buf;
    if (mpi_rank == 0) {
        ids_recv_buf = new int[total_num_scored_classifications];
        scores_recv_buf = new double[total_num_scored_classifications];
    }
    MPI_Gatherv(*ids,
                num_scored_classifications,
                MPI_INT,
                ids_recv_buf,
                mpi_counts,
                mpi_displs,
                MPI_INT,
                0,  // root
                MPI_COMM_WORLD);
    MPI_Gatherv(*scores,
                num_scored_classifications,
                MPI_DOUBLE,
                scores_recv_buf,
                mpi_counts,
                mpi_displs,
                MPI_DOUBLE,
                0,  // root
                MPI_COMM_WORLD);
    delete[] * ids;
    delete[] * scores;
    if (mpi_rank == 0) {
        *ids = ids_recv_buf;
        *scores = scores_recv_buf;
    }
}
#endif  // USE_MPI

void write_scored_classifications_raw(const string& filename,
                                      int num_scored_classifications,
                                      int* ids,
                                      double* scores) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    ofstream f(filename);
    f << "radmat_iz,best_pep_iz,best_pep_score\n";
    for (int i = 0; i < num_scored_classifications; i++) {
        f << i << ",";
        f << ids[i] << ",";
        f << setprecision(17) << scores[i] << "\n";
    }
    f.close();
    delete[] ids;
    delete[] scores;
}

}  // namespace fluoroseq
