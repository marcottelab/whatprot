// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "radiometries_io.h"

#include <fstream>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

#include "common/radiometry.h"
#ifdef USE_MPI
#include "io/mpi_counts_displs.h"
#endif  // USE_MPI

namespace fluoroseq {

namespace {
using std::ifstream;
using std::string;
}  // namespace

void read_radiometries(const string& filename,
                       int* num_timesteps,
                       int* num_channels,
                       int* total_num_radiometries,
                       int* num_radiometries,
                       Radiometry*** radiometries) {
    double* intensities;
    read_radiometries_raw(filename,
                          num_timesteps,
                          num_channels,
                          total_num_radiometries,
                          num_radiometries,
                          &intensities);
#ifdef USE_MPI
    scatter_radiometries(num_timesteps,
                         num_channels,
                         total_num_radiometries,
                         num_radiometries,
                         &intensities);
#else  // USE_MPI
    *num_radiometries = *total_num_radiometries;
#endif  // USE_MPI
    convert_radiometries_from_raw(*num_timesteps,
                                  *num_channels,
                                  *num_radiometries,
                                  intensities,
                                  radiometries);
    delete[] intensities;
}

void read_radiometries_raw(const string& filename,
                           int* num_timesteps,
                           int* num_channels,
                           int* total_num_radiometries,
                           int* num_radiometries,
                           double** intensities) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> *total_num_radiometries;
    *num_radiometries = *total_num_radiometries;
    *intensities = new double[(*total_num_radiometries)
                              * (*num_timesteps)
                              * (*num_channels)];
    for (   int i = 0;
            i < (*total_num_radiometries) * (*num_timesteps) * (*num_channels);
            i++) {
        f >> (*intensities)[i];
    }
    f.close();
}

#ifdef USE_MPI
void scatter_radiometries(int* num_timesteps,
                          int* num_channels,
                          int* total_num_radiometries,
                          int* num_radiometries,
                          double** intensities) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Bcast(num_timesteps,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(total_num_radiometries,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* mpi_counts;
    int* mpi_displs;
    mpi_counts_displs(*total_num_radiometries,
                      (*num_timesteps) * (*num_channels),  // block size
                      &mpi_counts,
                      &mpi_displs);
    int* mpi_num_radiometries;
    int* mpi_dummy_var;
    mpi_counts_displs(*total_num_radiometries,
                      1,  // block size
                      &mpi_num_radiometries,  // counts
                      &mpi_dummy_var);  // displs
    *num_radiometries = mpi_num_radiometries[mpi_rank];
    double* recv_buffer = new double[mpi_counts[mpi_rank]];
    MPI_Scatterv(*intensities,  // sendbuf
                 mpi_counts,  // sendcounts
                 mpi_displs,  // displs
                 MPI_DOUBLE,  // sendtype
                 recv_buffer,  // recvbuf
                 mpi_counts[mpi_rank],  // recvcount
                 MPI_DOUBLE,  // recvtype
                 0,  // root
                 MPI_COMM_WORLD);
    delete[] mpi_counts;
    delete[] mpi_displs;
    delete[] mpi_num_radiometries;
    delete[] mpi_dummy_var;
    if (mpi_rank == 0) {
        delete[] *intensities;
    }
    *intensities = recv_buffer;
}
#endif  // USE_MPI

void convert_radiometries_from_raw(int num_timesteps,
                                   int num_channels,
                                   int num_radiometries,
                                   double* intensities,
                                   Radiometry*** radiometries) {
    *radiometries = new Radiometry*[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        (*radiometries)[i] = new Radiometry(num_timesteps, num_channels);
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            (*radiometries)[i]->intensities[j] = intensities[
                    i * (num_timesteps * num_channels) + j];
        }
    }
}

}  // namespace fluoroseq
