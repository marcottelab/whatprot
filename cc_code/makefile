all: x64 skx
x64: classify classify_omp classify_mpi
skx: classify_skx classify_omp_skx classify_mpi_skx
classify: ./src/*/*
	g++ -O2 ./src/main/main.cc ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/tensor/*.cc -iquote ./src/ -o ./bin/classify
classify_omp: ./src/*/*
	g++ -fopenmp -O2 ./src/main/main.cc ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/tensor/*.cc -iquote ./src/ -o ./bin/classify_omp
classify_mpi: ./src/*/*
	mpic++ -O2 ./src/main/mpi_main.cc ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/tensor/*.cc -iquote ./src/ -o ./bin/classify_mpi
classify_skx: ./src/*/*
	icc -xCORE-AVX512 -O2 ./src/main/main.cc ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/tensor/*.cc -iquote ./src/ -o ./bin/classify_skx
classify_omp_skx: ./src/*/*
	icc -xCORE-AVX512 -qopenmp -O2 ./src/main/main.cc ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/tensor/*.cc -iquote ./src -o ./bin/classify_omp_skx
classify_mpi_skx: ./src/*/*
	mpicc -xCORE-AVX512 -O2 ./src/main/main.cc ./src/classifiers/*.cc ./src/common/*.cc ./src/fwd_alg/*.cc ./src/tensor/*.cc -iquote ./src/ -o ./bin/classify_mpi_skx
clean:
	rm -f ./bin/*
