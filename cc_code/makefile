# Author: Matthew Beauregard Smith (UT Austin)

OFLAG = -O2
CC_DEPS = ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/tensor/*.cc
IQUOTES = -iquote ./src/ -iquote ./extern/flann-1.9.1/src/cpp/
MAIN = ./src/main/main.cc
MPI_MAIN = ./src/main/mpi_main.cc

# LIBRARIES MUST BE AFTER NAMING ALL C++ FILES. BE CAREFUL WITH ORDER OF ARGUMENTS.
FLANN = -L ./extern/flann-1.9.1/build/lib/ -lflann_cpp_s

ARGS = $(OFLAG) $(MAIN) $(CC_DEPS) $(IQUOTES) $(FLANN)
MPI_ARGS = $(OFLAG) $(MPI_MAIN) $(CC_DEPS) $(IQUOTES) $(FLANN)

all: x64 skx
x64: classify classify_omp classify_mpi classify_mpi_omp
skx: classify_skx classify_omp_skx classify_mpi_skx classify_mpi_omp_skx
classify: ./src/*/*
	g++ $(ARGS) -o ./bin/classify
classify_omp: ./src/*/*
	g++ -fopenmp $(ARGS) -o ./bin/classify_omp
classify_mpi: ./src/*/*
	mpic++ $(MPI_ARGS) -o ./bin/classify_mpi
classify_mpi_omp: ./src/*/*
	mpic++ -fopenmp $(MPI_ARGS) -o ./bin/classify_mpi_omp
classify_skx: ./src/*/*
	icc -xCORE-AVX512 $(ARGS) -o ./bin/classify_skx
classify_omp_skx: ./src/*/*
	icc -xCORE-AVX512 -qopenmp $(ARGS) -o ./bin/classify_omp_skx
classify_mpi_skx: ./src/*/*
	mpicc -xCORE-AVX512 $(MPI_ARGS) -o ./bin/classify_mpi_skx
classify_mpi_omp_skx: ./src/*/*
	mpicc -xCORE-AVX512 -qopenmp $(MPI_ARGS) -o ./bin/classify_mpi_omp_skx
clean:
	rm -f ./bin/*
