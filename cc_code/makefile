# Author: Matthew Beauregard Smith (UT Austin)

OFLAG = -O2
# OFLAG = -Og -ggdb
CC_DEPS_SHARED = ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/io/*.cc ./src/main/*.cc ./src/simulation/*.cc ./src/tensor/*.cc ./src/util/*.cc
CC_DEPS_CLASSIFY = $(CC_DEPS_SHARED) ./src/main/classify/*.cc
CC_DEPS_SIMULATE = $(CC_DEPS_SHARED) ./src/main/simulate/*.cc
INCLUDES = -iquote ./src/ -I ./extern/flann-1.9.1/src/cpp/ -I ./extern/mrmpi-7Apr14/src/

# LIBRARIES MUST BE AFTER NAMING ALL C++ FILES. BE CAREFUL WITH ORDER OF ARGUMENTS.
FLANN = -L ./extern/flann-1.9.1/build/lib/ -lflann_cpp_s
MRMPI = -L ./extern/mrmpi-7Apr14/src/ -lmrmpi_mpicc
MRMPI_SERIAL = -L ./extern/mrmpi-7Apr14/src/ -lmrmpi_serial -I ./extern/mrmpi-7Apr14/mpistubs/ -L ./extern/mrmpi-7Apr14/mpistubs/ -lmpi_stubs
MRMPI_SKX = -L ./extern/mrmpi-7Apr14/src/ -lmrmpi_skx_mpicc
MRMPI_SKX_SERIAL = -L ./extern/mrmpi-7Apr14/src/ -lmrmpi_skx_serial -I ./extern/mrmpi-7Apr14/mpistubs/ -L ./extern/mrmpi-7Apr14/mpistubs/ -lmpi_stubs

ARGS_CLASSIFY = $(OFLAG) $(CC_DEPS_CLASSIFY) $(INCLUDES) $(FLANN) $(MRMPI_SERIAL)
MPI_ARGS_CLASSIFY = $(OFLAG) -D USE_MPI $(CC_DEPS_CLASSIFY) $(INCLUDES) $(FLANN) $(MRMPI)
ARGS_SIMULATE = $(OFLAG) $(CC_DEPS_SIMULATE) $(INCLUDES) $(FLANN) $(MRMPI_SERIAL)
MPI_ARGS_SIMULATE = $(OFLAG) -D USE_MPI $(CC_DEPS_SIMULATE) $(INCLUDES) $(FLANN) $(MRMPI)
ARGS_CLASSIFY_SKX = $(OFLAG) -wd3180 $(CC_DEPS_CLASSIFY) $(INCLUDES) $(FLANN) $(MRMPI_SKX_SERIAL)
MPI_ARGS_CLASSIFY_SKX = $(OFLAG) -wd3180 -D USE_MPI $(CC_DEPS_CLASSIFY) $(INCLUDES) $(FLANN) $(MRMPI_SKX)
ARGS_SIMULATE_SKX = $(OFLAG) -wd3180 $(CC_DEPS_SIMULATE) $(INCLUDES) $(FLANN) $(MRMPI_SKX_SERIAL)
MPI_ARGS_SIMULATE_SKX = $(OFLAG) -wd3180 -D USE_MPI $(CC_DEPS_SIMULATE) $(INCLUDES) $(FLANN) $(MRMPI_SKX)

all: x64 skx
x64: classify classify_mpi simulate simulate_mpi
skx: classify_skx classify_mpi_skx simulate_skx simulate_mpi_skx
classify: ./src/*/*
	g++ $(ARGS_CLASSIFY) -o ./bin/classify
classify_mpi: ./src/*/*
	mpic++ $(MPI_ARGS_CLASSIFY) -o ./bin/classify_mpi
classify_skx: ./src/*/*
	icc -xCORE-AVX512 $(ARGS_CLASSIFY_SKX) -o ./bin/classify_skx
classify_mpi_skx: ./src/*/*
	mpicc -xCORE-AVX512 $(MPI_ARGS_CLASSIFY_SKX) -o ./bin/classify_mpi_skx
simulate: ./src/*/*
	g++ $(ARGS_SIMULATE) -o ./bin/simulate
simulate_mpi: ./src/*/*
	mpic++ $(MPI_ARGS_SIMULATE) -o ./bin/simulate_mpi
simulate_skx: ./src/*/*
	icc -xCORE-AVX512 $(ARGS_SIMULATE_SKX) -o ./bin/simulate_skx
simulate_mpi_skx: ./src/*/*
	mpicc -xCORE-AVX512 $(MPI_ARGS_SIMULATE_SKX) -o ./bin/simulate_mpi_skx
clean:
	rm -f ./bin/*
