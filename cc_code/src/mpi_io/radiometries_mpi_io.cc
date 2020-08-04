// Author: Matthew Beauregard Smith (UT Austin)
#include "radiometries_mpi_io.h"

#include <fstream>
#include <string>

#include <mpi.h>

#include "common/radiometry.h"

namespace fluoroseq {

namespace {
using std::ifstream;
using std::string;
}  // namespace

void mpi_read_radiometries(const std::string& filename,
                           int* num_timesteps,
                           int* num_channels,
                           int* total_num_radiometries,
                           int* num_radiometries,
                           Radiometry*** radiometries) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == 0) {
        mpi_read_radiometries_master(filename,
                                     num_timesteps,
                                     num_channels,
                                     total_num_radiometries,
                                     num_radiometries,
                                     radiometries);
    } else {
        mpi_read_radiometries_slave(filename,
                                    num_timesteps,
                                    num_channels,
                                    total_num_radiometries,
                                    num_radiometries,
                                    radiometries);
    }
}

void mpi_read_radiometries_master(const std::string& filename,
                                  int* num_timesteps,
                                  int* num_channels,
                                  int* total_num_radiometries,
                                  int* num_radiometries,
                                  Radiometry*** radiometries) {
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    ifstream f(filename);
    f >> *num_timesteps;
    MPI_Bcast(num_timesteps,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    f >> *num_channels;
    MPI_Bcast(num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int total_num_radiometries;
    f >> *total_num_radiometries;
    MPI_Bcast(&total_num_radiometries,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    double* intensities = new double[(*total_num_radiometries)
                                     * (*num_timesteps)
                                     * (*num_channels)];
    for (   int i = 0;
            i < (*total_num_radiometries) * (*num_timesteps) * (*num_channels);
            i++) {
        f >> intensities[i];
    }
    int* intensity_sendcounts = new int[mpi_size];
    int* intensity_displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long) (*total_num_radiometries) * (long) i
                    / (long) mpi_size;
        int end = (long) (*total_num_radiometries) * (long) (i + 1)
                  / (long) mpi_size;
        intensity_sendcounts[i] = (end - begin)
                                  * (*num_timesteps)
                                  * (*num_channels);
        intensity_displs[i] = begin
                              * (*num_timesteps)
                              * (*num_channels);
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
    (*num_radiometries) = intensity_sendcounts[0]
                          / (*num_timesteps)
                          / (*num_channels);
    *radiometries = new Radiometry*[*num_radiometries];
    for (int i = 0; i < *num_radiometries; i++) {
        (*radiometries)[i] = new Radiometry(*num_timesteps, *num_channels);
        for (int j = 0; j < *num_timesteps * (*num_channels); j++) {
            (*radiometries)[i]->intensities[j] = intensities[
                    i * ((*num_timesteps) * (*num_channels)) + j];
        }
    }
    f.close();
    delete[] intensities;
    delete[] intensity_sendcounts;
    delete[] intensity_displs;
    delete[] dummy_recv_buffer;
}

void mpi_read_radiometries_slave(const std::string& filename,
                                 int* num_timesteps,
                                 int* num_channels,
                                 int* total_num_radiometries,
                                 int* num_radiometries,
                                 Radiometry*** radiometries) {
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
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
    int* intensity_sendcounts = new int[mpi_size];
    int* intensity_displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long) (*total_num_radiometries) * (long) i / (long) mpi_size;
        int end = (long) (*total_num_radiometries) * (long) (i + 1)
                  / (long) mpi_size;
        intensity_sendcounts[i] = (end - begin) * (*num_timesteps) * (*num_channels);
        intensity_displs[i] = begin * (*num_timesteps) * (*num_channels);
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
                           / (*num_timesteps)
                           / (*num_channels);
    Radiometry** radiometries = new Radiometry*[*num_radiometries];
    for (int i = 0; i < *num_radiometries; i++) {
        (*radiometries)[i] = new Radiometry(*num_timesteps, *num_channels);
        for (int j = 0; j < (*num_timesteps) * (*num_channels); j++) {
            (*radiometries)[i]->intensities[j] = intensities[
                    i * ((*num_timesteps) * (*num_channels)) + j];
        }
    }
    delete[] intensity_sendcounts;
    delete[] intensity_displs;
    delete[] intensities;
}

}  // namespace fluoroseq
