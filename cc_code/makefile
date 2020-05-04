all: classify classify_omp
classify: ./src/*/*
	g++ -O2 ./src/*/*.cc -iquote ./src/ -o ./bin/classify
classify_omp: ./src/*/*
	g++ -fopenmp -O2 ./src/*/*.cc -iquote ./src/ -o ./bin/classify_omp
clean:
	rm ./bin/classify
	rm ./bin/classify_omp