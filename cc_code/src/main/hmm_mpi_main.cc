// Author: Matthew Beauregard Smith (UT Austin)
// Simple application to read in a TSV file of dye seqs, a TSV file of
// radiometries, and write predicted classifications to a TSV file.

#include <cstddef>  // for offsetof
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <iostream>
#include <string>

// Not a standard c++ library. Must have installed MPI to build.
#include <mpi.h>

#include "common/approximation_model.h"
#include "common/dye_seq.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "util/time.h"

#ifdef _OPENMP
#include "classifiers/omp_fwd_alg_classifier.h"
#else
#include "classifiers/fwd_alg_classifier.h"
#endif

namespace {
using fluoroseq::wall_time;
using fluoroseq::ApproximationModel;
using fluoroseq::DistributionType;
using fluoroseq::DyeSeq;
using fluoroseq::ErrorModel;
using fluoroseq::Radiometry;
using fluoroseq::ScoredClassification;
using fluoroseq::SourcedData;
using fluoroseq::SourceCount;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::setprecision;
using std::string;
#ifdef _OPENMP
using fluoroseq::OMPFwdAlgClassifier;
#else
using fluoroseq::FwdAlgClassifier;
#endif
}

void main_for_master(char** argv);
void main_for_slave();

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (argc != 4) {
        if (mpi_rank == 0) {
            cout << "wrong number of inputs\n";
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (mpi_rank == 0) {
        main_for_master(argv);
    } else {
        main_for_slave();
    }
    MPI_Finalize();
    return 0;
}

void main_for_master(char** argv) {
    double total_start_time = wall_time();

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    cout << "Using MPI\n";
    cout << "    Number of processes: " << mpi_size << "\n";

    #ifdef _OPENMP
    cout << "Using OpenMP\n";
    cout << "    Threads per process: " << omp_get_max_threads() << "\n";
    #endif

    char* dye_seqs_filename = argv[1];
    char* radiometries_filename = argv[2];
    char* predictions_filename = argv[3];

    double start_time = wall_time();
    ErrorModel error_model(.06,  // p_edman_failure
                           .05,  // p_detach
                           .05,  // p_bleach
                           .07,  // p_dud
                           DistributionType::LOGNORMAL,
                           1.0,  // mu
                           .16);  // sigma
    ApproximationModel approximation_model(16);
    double end_time = wall_time();
    cout << "Built error model.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    start_time = wall_time();
    ifstream fdye(dye_seqs_filename);
    int num_channels;
    fdye >> num_channels;
    MPI_Bcast(&num_channels,
               1,  // count
               MPI_INT,
               0,  // root
               MPI_COMM_WORLD);
    int num_dye_seqs;
    fdye >> num_dye_seqs;
    MPI_Bcast(&num_dye_seqs,
               1, // count
               MPI_INT,
               0,  // root
               MPI_COMM_WORLD);
    string* dye_strings = new string[num_dye_seqs];
    int* dye_string_lengths = new int[num_dye_seqs];
    int* dye_seqs_num_peptides = new int[num_dye_seqs];
    int* dye_seqs_ids = new int[num_dye_seqs];
    for (int i = 0; i < num_dye_seqs; i++) {
        fdye >> dye_strings[i];
        dye_string_lengths[i] = dye_strings[i].length();
        fdye >> dye_seqs_num_peptides[i];
        fdye >> dye_seqs_ids[i];
    }
    fdye.close();
    MPI_Bcast(dye_string_lengths,
              num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    for (int i = 0; i < num_dye_seqs; i++) {
        MPI_Bcast(const_cast<char *>(dye_strings[i].c_str()),
                  dye_string_lengths[i],
                  MPI_CHAR,
                  0,  // root
                  MPI_COMM_WORLD);
    }
    MPI_Bcast(dye_seqs_num_peptides,
              num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(dye_seqs_ids,
              num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    SourcedData<DyeSeq*, SourceCount<int>>** dye_seqs; 
    dye_seqs = new SourcedData<DyeSeq*, SourceCount<int>>*[num_dye_seqs];
    for (int i = 0; i < num_dye_seqs; i++) {
        dye_seqs[i] = new SourcedData<DyeSeq*, SourceCount<int>>(
                              new DyeSeq(num_channels, dye_strings[i]),
                              SourceCount<int>(dye_seqs_ids[i],
                                                   dye_seqs_num_peptides[i]));
    }
    end_time = wall_time();
    cout << "Read from dye seqs file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of dye seqs: " << num_dye_seqs << "\n";

    start_time = wall_time();
    ifstream frad(radiometries_filename);
    int num_timesteps;
    frad >> num_timesteps;
    MPI_Bcast(&num_timesteps,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    // Note, already have num channels. Should be writing over this with same
    // value. Therefore don't need to bcast it with MPI.
    int duplicate_num_channels;
    frad >> duplicate_num_channels;
    int total_num_radiometries;
    frad >> total_num_radiometries;
    MPI_Bcast(&total_num_radiometries,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    double* intensities = new double[total_num_radiometries
                                     * num_timesteps
                                     * num_channels];
    for (   int i = 0;
            i < total_num_radiometries * num_timesteps * num_channels;
            i++) {
        frad >> intensities[i];
    }
    int* intensity_sendcounts = new int[mpi_size];
    int* intensity_displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long) total_num_radiometries * (long) i / (long) mpi_size;
        int end = (long) total_num_radiometries * (long) (i + 1)
                  / (long) mpi_size;
        intensity_sendcounts[i] = (end - begin) * num_timesteps * num_channels;
        intensity_displs[i] = begin * num_timesteps * num_channels;
    }
    double* dummy_recv_buffer = new double[intensity_sendcounts[0]];
    MPI_Scatterv(intensities,
                 intensity_sendcounts,
                 intensity_displs,
                 MPI_DOUBLE,
                 dummy_recv_buffer,
                 intensity_sendcounts[0],
                 MPI_DOUBLE,
                 0,  // root
                 MPI_COMM_WORLD);
    int num_radiometries = intensity_sendcounts[0]
                           / num_timesteps
                           / num_channels;
    Radiometry** radiometries = new Radiometry*[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        radiometries[i] = new Radiometry(num_timesteps, num_channels);
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            radiometries[i]->intensities[j] = intensities[i * (num_timesteps
                                                               * num_channels)
                                                          + j];
        }
    }
    frad.close();
    end_time = wall_time();
    cout << "Read from radiometries file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Number of radiometries: " << num_radiometries << "\n";

    start_time = wall_time();
    #ifdef _OPENMP
    OMPFwdAlgClassifier classifier(num_timesteps,
                                   num_channels,
                                   error_model,
                                   approximation_model,
                                   num_dye_seqs,
                                   dye_seqs);
    #else
    FwdAlgClassifier classifier(num_timesteps,
                                num_channels,
                                error_model,
                                approximation_model,
                                num_dye_seqs,
                                dye_seqs);
    #endif
    end_time = wall_time();
    cout << "Constructed classifier.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    start_time = wall_time();
    ScoredClassification* results = classifier.classify(num_radiometries,
                                                        radiometries);
    ScoredClassification* total_results = new ScoredClassification[
            total_num_radiometries];
    int* result_sendcounts = new int[mpi_size];
    int* result_displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long) total_num_radiometries * (long) i / (long) mpi_size;
        int end = (long) total_num_radiometries * (long) (i + 1)
                  / (long) mpi_size;
        result_sendcounts[i] = end - begin;
        result_displs[i] = begin;
    }
    int array_of_blocklengths[] = {1, 1, 1};
    MPI_Aint array_of_displacements[] = {
            (MPI_Aint) offsetof(ScoredClassification, score),
            (MPI_Aint) offsetof(ScoredClassification, total),
            (MPI_Aint) offsetof(ScoredClassification, id)};
    MPI_Datatype array_of_types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Datatype MPI_SCORED_CLASSIFICATION_NOT_RESIZED;
    MPI_Type_struct(3,  // count
                    array_of_blocklengths,
                    array_of_displacements,
                    array_of_types,
                    &MPI_SCORED_CLASSIFICATION_NOT_RESIZED);
    MPI_Aint lb;
    MPI_Aint extent;
    MPI_Type_get_extent(MPI_SCORED_CLASSIFICATION_NOT_RESIZED, &lb, &extent);
    MPI_Datatype MPI_SCORED_CLASSIFICATION;
    MPI_Type_create_resized(MPI_SCORED_CLASSIFICATION_NOT_RESIZED,
                            lb,
                            extent,
                            &MPI_SCORED_CLASSIFICATION);
    MPI_Type_commit(&MPI_SCORED_CLASSIFICATION);
    MPI_Type_free(&MPI_SCORED_CLASSIFICATION_NOT_RESIZED);
    MPI_Gatherv(results,
                num_radiometries,
                MPI_SCORED_CLASSIFICATION,
                total_results,
                result_sendcounts,
                result_displs,
                MPI_SCORED_CLASSIFICATION,
                0,  // root
                MPI_COMM_WORLD);
    MPI_Type_free(&MPI_SCORED_CLASSIFICATION);
    end_time = wall_time();
    cout << "Got results.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";
    cout << "    Total number of runs: "
        << num_dye_seqs * num_radiometries
        << "\n";
    cout << "    Time per run in seconds: "
        << (end_time - start_time) / (num_dye_seqs * num_radiometries)
        << "\n";

    start_time = wall_time();
    ofstream fpred(predictions_filename);
    fpred << "radmat_iz,best_pep_iz,best_pep_score\n";
    for (int i = 0; i < total_num_radiometries; i++) {
        fpred << i << ",";
        fpred << total_results[i].id << ",";
        fpred << setprecision(17) << total_results[i].adjusted_score() << "\n";
    }
    fpred.close();
    end_time = wall_time();
    cout << "Wrote to file.\n";
    cout << "    Time in seconds: " << end_time - start_time << "\n";

    double total_end_time = wall_time();
    cout << "Totals:\n";
    cout << "    Time in seconds: "
        << total_end_time - total_start_time
        << "\n";
    cout << "    Total number of runs: "
        << num_dye_seqs * num_radiometries
        << "\n";
    cout << "    Time per run in seconds: "
        << (total_end_time - total_start_time)
            / (num_dye_seqs * num_radiometries)
        << "\n";
}

void main_for_slave() {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ErrorModel error_model(.06,  // p_edman_failure
                           .05,  // p_detach
                           .05,  // p_bleach
                           .07,  // p_dud
                           DistributionType::LOGNORMAL,
                           1.0,  // mu
                           .16);  // sigma
    ApproximationModel approximation_model(16);
    
    int num_channels;
    MPI_Bcast(&num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int num_dye_seqs;
    MPI_Bcast(&num_dye_seqs,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* dye_string_lengths = new int[num_dye_seqs];
    MPI_Bcast(dye_string_lengths,
              num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    string* dye_strings = new string[num_dye_seqs];
    for (int i = 0; i < num_dye_seqs; i++) {
        char* cstr = new char[dye_string_lengths[i]];
        MPI_Bcast(cstr,
                  dye_string_lengths[i],
                  MPI_CHAR,
                  0,  // root
                  MPI_COMM_WORLD);
        dye_strings[i] = string(cstr, dye_string_lengths[i]);
        delete[] cstr;
    }
    delete[] dye_string_lengths;
    int* dye_seqs_num_peptides = new int[num_dye_seqs];
    MPI_Bcast(dye_seqs_num_peptides,
              num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* dye_seqs_ids = new int[num_dye_seqs];
    MPI_Bcast(dye_seqs_ids,
              num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    SourcedData<DyeSeq*, SourceCount<int>>** dye_seqs; 
    dye_seqs = new SourcedData<DyeSeq*, SourceCount<int>>*[num_dye_seqs];
    for (int i = 0; i < num_dye_seqs; i++) {
        dye_seqs[i] = new SourcedData<DyeSeq*, SourceCount<int>>(
                              new DyeSeq(num_channels, dye_strings[i]),
                              SourceCount<int>(dye_seqs_ids[i],
                                                   dye_seqs_num_peptides[i]));
    }

    int num_timesteps;
    MPI_Bcast(&num_timesteps,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int total_num_radiometries;
    MPI_Bcast(&total_num_radiometries,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* intensity_sendcounts = new int[mpi_size];
    int* intensity_displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long) total_num_radiometries * (long) i / (long) mpi_size;
        int end = (long) total_num_radiometries * (long) (i + 1)
                  / (long) mpi_size;
        intensity_sendcounts[i] = (end - begin) * num_timesteps * num_channels;
        intensity_displs[i] = begin * num_timesteps * num_channels;
    }
    double* intensities = new double[intensity_sendcounts[mpi_rank]];
    MPI_Scatterv(NULL,  // sendbuf
                 intensity_sendcounts,
                 intensity_displs,
                 MPI_DOUBLE,
                 intensities,
                 intensity_sendcounts[mpi_rank],
                 MPI_DOUBLE,
                 0,  // root
                 MPI_COMM_WORLD);
    int num_radiometries = intensity_sendcounts[mpi_rank]
                           / num_timesteps
                           / num_channels;
    Radiometry** radiometries = new Radiometry*[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        radiometries[i] = new Radiometry(num_timesteps, num_channels);
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            radiometries[i]->intensities[j] = intensities[i * (num_timesteps
                                                               * num_channels)
                                                          + j];
        }
    }

    #ifdef _OPENMP
    OMPFwdAlgClassifier classifier(num_timesteps,
                                   num_channels,
                                   error_model,
                                   approximation_model,
                                   num_dye_seqs,
                                   dye_seqs);
    #else
    FwdAlgClassifier classifier(num_timesteps,
                                num_channels,
                                error_model,
                                approximation_model,
                                num_dye_seqs,
                                dye_seqs);
    #endif

    ScoredClassification* results = classifier.classify(num_radiometries,
                                                        radiometries);
    
    int* result_sendcounts = new int[mpi_size];
    int* result_displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long) total_num_radiometries * (long) i / (long) mpi_size;
        int end = (long) total_num_radiometries * (long) (i + 1)
                  / (long) mpi_size;
        result_sendcounts[i] = end - begin;
        result_displs[i] = begin;
    }
    int array_of_blocklengths[] = {1, 1, 1};
    MPI_Aint array_of_displacements[] = {
            (MPI_Aint) offsetof(ScoredClassification, score),
            (MPI_Aint) offsetof(ScoredClassification, total),
            (MPI_Aint) offsetof(ScoredClassification, id)};
    MPI_Datatype array_of_types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Datatype MPI_SCORED_CLASSIFICATION_NOT_RESIZED;
    MPI_Type_struct(3,  // count
                    array_of_blocklengths,
                    array_of_displacements,
                    array_of_types,
                    &MPI_SCORED_CLASSIFICATION_NOT_RESIZED);
    MPI_Aint lb;
    MPI_Aint extent;
    MPI_Type_get_extent(MPI_SCORED_CLASSIFICATION_NOT_RESIZED, &lb, &extent);
    MPI_Datatype MPI_SCORED_CLASSIFICATION;
    MPI_Type_create_resized(MPI_SCORED_CLASSIFICATION_NOT_RESIZED,
                            lb,
                            extent,
                            &MPI_SCORED_CLASSIFICATION);
    MPI_Type_commit(&MPI_SCORED_CLASSIFICATION);
    MPI_Type_free(&MPI_SCORED_CLASSIFICATION_NOT_RESIZED);
    MPI_Gatherv(results,
                num_radiometries,
                MPI_SCORED_CLASSIFICATION,
                NULL,  // recvbuf
                NULL,  // recvcounts
                NULL,  // displs
                MPI_SCORED_CLASSIFICATION,
                0,  // root
                MPI_COMM_WORLD);
    MPI_Type_free(&MPI_SCORED_CLASSIFICATION);
}
