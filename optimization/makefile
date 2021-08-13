.PHONY: all
all: cell_distances

cell_distances: cell_distances.c
	gcc -O3 -fopenmp -o cell_distances cell_distances.c -lm -lgomp

omp_test: omp_test.c
	gcc -O2 -fopenmp -o omp_test omp_test.c -lm -lgomp

run: cell_distances
	./cell_distances -t5

cell_distances.tar.gz: cell_distances.c makefile
	tar -cvzf cell_distances.tar.gz cell_distances.c makefile

.PHONY: test
test: clean cell_distances.tar.gz
	/home/hpc2020/cell_distances/check_submission.py /home/hpcuser319/HPC_ASSIGN2/cell_distances.tar.gz

.PHONY: clean
clean:
	rm -rf cell_distances distances/ extracted/ distances.tar.gz
