/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// For MPI version, define compiler macro USE_MPI when building.

// Defining symbols from header:
#include "radiometries-io.h"

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
#include "common/radiometry.h"
#ifdef USE_MPI
#include "io/mpi-counts-displs.h"
#endif  // USE_MPI

namespace fluoroseq {

namespace {
using std::ifstream;
using std::ofstream;
using std::setprecision;
using std::string;
using std::vector;
}  // namespace

void read_radiometries(const string& filename,
                       int* num_timesteps,
                       int* num_channels,
                       int* total_num_radiometries,
                       vector<Radiometry>* radiometries) {
    int num_radiometries;
    double* intensities;
    read_radiometries_raw(filename,
                          num_timesteps,
                          num_channels,
                          total_num_radiometries,
                          &num_radiometries,
                          &intensities);
#ifdef USE_MPI
    scatter_radiometries(num_timesteps,
                         num_channels,
                         total_num_radiometries,
                         &num_radiometries,
                         &intensities);
#else  // USE_MPI
    num_radiometries = *total_num_radiometries;
#endif  // USE_MPI
    convert_radiometries_from_raw(*num_timesteps,
                                  *num_channels,
                                  num_radiometries,
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
    *intensities = new double[(*total_num_radiometries) * (*num_timesteps)
                              * (*num_channels)];
    for (int i = 0;
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
        delete[] * intensities;
    }
    *intensities = recv_buffer;
}
#endif  // USE_MPI

void convert_radiometries_from_raw(int num_timesteps,
                                   int num_channels,
                                   int num_radiometries,
                                   double* intensities,
                                   vector<Radiometry>* radiometries) {
    radiometries->reserve(num_radiometries);
    for (int i = 0; i < num_radiometries; i++) {
        radiometries->push_back(Radiometry(num_timesteps, num_channels));
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            radiometries->back().intensities[j] =
                    intensities[i * (num_timesteps * num_channels) + j];
        }
    }
}

void write_radiometries(
        const string& filename,
        int num_timesteps,
        int num_channels,
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries) {
    double* intensities;
    convert_raw_from_radiometries(
            radiometries,
            num_timesteps * num_channels,  // radiometry size
            &intensities);
#ifdef USE_MPI
    int total_num_radiometries;
    gather_radiometries(radiometries.size(),  // number of radiometries
                        num_timesteps * num_channels,  // radiometry size
                        &intensities,
                        &total_num_radiometries);
#else  // USE_MPI
    int total_num_radiometries = radiometries.size();
#endif  // USE_MPI
    write_radiometries_raw(filename,
                           num_timesteps,
                           num_channels,
                           total_num_radiometries,
                           intensities);
}

void convert_raw_from_radiometries(
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries,
        int radiometry_size,
        double** intensities) {
    *intensities = new double[radiometries.size() * radiometry_size];
    for (int i = 0; i < radiometries.size(); i++) {
        for (int j = 0; j < radiometry_size; j++) {
            (*intensities)[i * radiometry_size + j] =
                    radiometries[i].value.intensities[j];
        }
    }
}

#ifdef USE_MPI
void gather_radiometries(int num_radiometries,
                         int radiometry_size,
                         double** intensities,
                         int* total_num_radiometries) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    int* mpi_counts;
    if (mpi_rank == 0) {
        mpi_counts = new int[mpi_size];
    }
    MPI_Gather(&num_radiometries,  // sendbuf
               1,  // sendcount
               MPI_INT,
               mpi_counts,  // recvbuf
               1,  // recvcount, this is PER PROCESS
               MPI_INT,
               0,  // root
               MPI_COMM_WORLD);
    int* mpi_displs;
    if (mpi_rank == 0) {
        mpi_displs = new int[mpi_size];
        int tally = 0;
        for (int i = 0; i < mpi_size; i++) {
            mpi_displs[i] = tally * radiometry_size;
            tally += mpi_counts[i];
            mpi_counts[i] *= radiometry_size;
        }
        *total_num_radiometries = tally;
    }
    double* intensities_recv_buf;
    if (mpi_rank == 0) {
        intensities_recv_buf =
                new double[(*total_num_radiometries) * radiometry_size];
    }
    MPI_Gatherv(*intensities,  // sendbuf
                num_radiometries * radiometry_size,  // sendcount
                MPI_DOUBLE,
                intensities_recv_buf,  // recvbuf
                mpi_counts,
                mpi_displs,
                MPI_DOUBLE,
                0,  // root
                MPI_COMM_WORLD);
    delete[] * intensities;
    if (mpi_rank == 0) {
        *intensities = intensities_recv_buf;
        delete[] mpi_counts;
        delete[] mpi_displs;
    }
}
#endif  // USE_MPI

void write_radiometries_raw(const std::string& filename,
                            int num_timesteps,
                            int num_channels,
                            int num_radiometries,
                            double* intensities) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    ofstream f(filename);
    f << num_timesteps << "\n";
    f << num_channels << "\n";
    f << num_radiometries << "\n";
    for (int i = 0; i < num_radiometries; i++) {
        for (int j = 0; j < num_timesteps * num_channels; j++) {
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
#ifdef USE_MPI
    int total_num_radiometries;
    gather_ys(radiometries.size(), &ys, &total_num_radiometries);
#else  // USE_MPI
    int total_num_radiometries = radiometries.size();
#endif  // USE_MPI
    write_ys_raw(filename,
                 total_num_radiometries,  // num radiometries
                 ys);
}

void get_raw_ys(
        const vector<SourcedData<Radiometry, SourceCount<int>>>& radiometries,
        int** ys) {
    *ys = new int[radiometries.size()];
    for (int i = 0; i < radiometries.size(); i++) {
        (*ys)[i] = radiometries[i].source.source;
    }
}

#ifdef USE_MPI
void gather_ys(int num_radiometries, int** ys, int* total_num_radiometries) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    int* mpi_counts;
    if (mpi_rank == 0) {
        mpi_counts = new int[mpi_size];
    }
    MPI_Gather(&num_radiometries,  // sendbuf
               1,  // sendcount
               MPI_INT,
               mpi_counts,  // recvbuf
               1,  // recvcount, this is PER PROCESS
               MPI_INT,
               0,  // root
               MPI_COMM_WORLD);
    int* mpi_displs;
    if (mpi_rank == 0) {
        mpi_displs = new int[mpi_size];
        int tally = 0;
        for (int i = 0; i < mpi_size; i++) {
            mpi_displs[i] = tally;
            tally += mpi_counts[i];
        }
        *total_num_radiometries = tally;
    }
    int* ys_recv_buf;
    if (mpi_rank == 0) {
        ys_recv_buf = new int[*total_num_radiometries];
    }
    MPI_Gatherv(*ys,
                num_radiometries,  // sendcount
                MPI_INT,
                ys_recv_buf,
                mpi_counts,
                mpi_displs,
                MPI_INT,
                0,  // root
                MPI_COMM_WORLD);
    delete[] * ys;
    if (mpi_rank == 0) {
        *ys = ys_recv_buf;
        delete[] mpi_counts;
        delete[] mpi_displs;
    }
}
#endif  // USE_MPI

void write_ys_raw(const string& filename, int num_radiometries, int* ys) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    ofstream f(filename);
    f << num_radiometries << "\n";
    for (int i = 0; i < num_radiometries; i++) {
        f << ys[i] << "\n";
    }
    f.close();
    delete[] ys;
}

}  // namespace fluoroseq
