# Author: Matthew Beauregard Smith (UT Austin)

OFLAG = -O2
# OFLAG = -Og -ggdb
CC_DEPS = ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/io/*.cc ./src/main/*.cc ./src/main/classify/*.cc ./src/tensor/*.cc ./src/util/*.cc
INCLUDES = -iquote ./src/ -I ./extern/flann-1.9.1/src/cpp/

# LIBRARIES MUST BE AFTER NAMING ALL C++ FILES. BE CAREFUL WITH ORDER OF ARGUMENTS.
FLANN = -L ./extern/flann-1.9.1/build/lib/ -lflann_cpp_s

ARGS = $(OFLAG) $(CC_DEPS) $(INCLUDES) $(FLANN)
MPI_ARGS = $(OFLAG) -D USE_MPI $(CC_DEPS) $(INCLUDES) $(FLANN)

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
