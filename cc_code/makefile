# Author: Matthew Beauregard Smith (UT Austin)

OFLAG = -O2
# OFLAG = -Og -ggdb
CC_DEPS = ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/io/*.cc ./src/tensor/*.cc ./src/util/range.cc
INCLUDES = -iquote ./src/ -I ./extern/flann-1.9.1/src/cpp/
MAIN = ./src/main/classify_main.cc ./src/main/ann_main.cc ./src/main/hmm_main.cc ./src/main/hybrid_main.cc ./src/util/time.cc
MPI_MAIN = ./src/main/hmm_mpi_main.cc ./src/util/mpi_time.cc

# LIBRARIES MUST BE AFTER NAMING ALL C++ FILES. BE CAREFUL WITH ORDER OF ARGUMENTS.
FLANN = -L ./extern/flann-1.9.1/build/lib/ -lflann_cpp_s

ARGS = $(OFLAG) $(MAIN) $(CC_DEPS) $(INCLUDES) $(FLANN)
MPI_ARGS = $(OFLAG) $(MPI_MAIN) $(CC_DEPS) $(INCLUDES) $(FLANN)

all: x64 skx
x64: classify classify_mpi
skx: classify_skx classify_mpi_skx
classify: ./src/*/*
	g++ $(ARGS) -o ./bin/classify
classify_mpi: ./src/*/*
	mpic++ $(MPI_ARGS) -o ./bin/classify_mpi
classify_skx: ./src/*/*
	icc -xCORE-AVX512 $(ARGS) -o ./bin/classify_skx
classify_mpi_skx: ./src/*/*
	mpicc -xCORE-AVX512 $(MPI_ARGS) -o ./bin/classify_mpi_skx
clean:
	rm -f ./bin/*
