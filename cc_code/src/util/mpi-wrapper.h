/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_UTIL_MPI_WRAPPER_H
#define WHATPROT_UTIL_MPI_WRAPPER_H

// Standard C++ library headers:
#include <sstream>
#include <vector>

// MPI
#include <mpi.h>

// External headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

namespace whatprot {

int mpi_size() {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    return world_size;
}

int mpi_rank() {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return world_rank;
}

template<class T>
void scatterv(const std::vector<T>& in_values, std::vector<T>* out_values, int root) {
    if (mpi_rank() == root) {
        int displ = 0;
        for (int dest = 0; dest < mpi_size(); dest++) {
            int count = (in_values.size() - displ) / (mpi_size() - dest);
            if (dest == root) {
                for (int i = 0; i < count; i++) {
                    out_values->push_back(in_values[displ + i]);
                }
            } else {
                std::stringstream ss;
                {  // unnamed scope to force oa to flush to ss.
                    cereal::BinaryOutputArchive oa(ss);
                    std::vector<T> dest_out_values;
                    for (int i = 0; i < count; i++) {
                        dest_out_values.push_back(in_values[displ + i]);
                    }
                    oa(dest_out_values);
                }
                std::string str = ss.str();
                MPI_Send(str.c_str(), str.length(), MPI_CHAR, dest, 0, MPI_COMM_WORLD);
            }
            displ += count;
        }
    } else {
        MPI_Status status;
        MPI_Probe(root, 0, MPI_COMM_WORLD, &status);
        int count;
        MPI_Get_count(&status, MPI_CHAR, &count);
        char* buf = new char[count];
        MPI_Recv(buf, count, MPI_CHAR, root, 0, MPI_COMM_WORLD);
        std::string str(buf, count);
        delete[] buf;
        std::stringstream ss(str);
        {  // not clear if necessary but online example for cereal scopes this.
            cereal::BinaryInputArchive ia(ss);
            ia(*out_values);
        }
    }
}

template<class T>
void gatherv(const std::vector<T>& in_values, std::vector<T>* out_values, int root) {
    if (mpi_rank() == root) {
        
    }
}

}  // namespace whatprot

#endif WHATPROT_UTIL_MPI_WRAPPER_H