all: x64 skx
x64: classify classify_omp
skx: classify_skx classify_omp_skx
classify: ./src/*/*
	g++ -O2 ./src/*/*.cc -iquote ./src/ -o ./bin/classify
classify_omp: ./src/*/*
	g++ -fopenmp -O2 ./src/*/*.cc -iquote ./src/ -o ./bin/classify_omp
classify_skx: ./src/*/*
	icc -xCORE-AVX512 -O2 ./src/*/*.cc -iquote ./src/ -o ./bin/classify_skx
classify_omp_skx: ./src/*/*
	icc -xCORE-AVX512 -qopenmp -O2 ./src/*/*.cc -iquote ./src -o ./bin/classify_omp_skx
clean:
	rm -f ./bin/*
