.PHONY: all
all: newton

newton: newton.c
	gcc -o newton newton.c -O2 -lpthread -lm
.PHONY: images
images: newton
	for d in {0..9}; do \
	echo "d=$$d";\
	./newton -t5 -l1000 $$d;\
	done
newton.tar.gz: newton.c makefile
	tar -cvzf newton.tar.gz newton.c makefile

.PHONY: test
test: clean newton.tar.gz
	/home/hpc2020/newton_iteration/check_submission.py /home/hpcuser319/hpc_assign3/newton.tar.gz


.PHONY: clean
clean:
	rm -rf newton extracted/ newton.tar.gz
