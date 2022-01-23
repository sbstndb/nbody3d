
all: aos1.g aos2.g soa1.g soa2.g soa3.g soa4.g soa5.g soatest.g

aos1.g:src/aos1.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp
	
	
aos2.g: src/aos2.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -Ofast $< -o $@ -lm -fopenmp
	

soa1.g:src/soa1.c
	gcc -march=native -mtune=native -fno-unroll-loops -fno-inline-functions -fno-peel-loops -fno-tree-vectorize -fno-tree-loop-vectorize -Ofast $< -o $@ -lm -fopenmp
	
soa2.g:src/soa2.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp
	
soa3.g:src/soa3.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp
	
soa4.g:src/soa4.c
	gcc -march=native -mtune=native -funroll-loops -funroll-all-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -mfma -fno-math-errno -Ofast $< -o $@ -lm -fopenmp			

soa5.g:src/soa5.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp


soatest.g:src/soatest.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp


aos1.i: src/aos1.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

clean:
	rm -Rf *~ *.g *.i *.optrpt *.txt
	
run:
	./aos1.g && ./aos2.g &&  ./soa1.g && ./soa2.g && ./soa3.g && ./soa4.g && ./soa5.g
	
	
run_one_512:
	 ./aos1.g 16000 0 1 && ./aos2.g 16000 0 1  &&  ./soa1.g 16000 0 1  && ./soa2.g 16000 0 1  && ./soa3.g 16000 0 1  && ./soa4.g 16000 0 1  && ./soa5.g 16000 0 1   && python3 error512.py >> result_precision_one.txt && cat result_precision_one.txt
	
run_one_2:
	./aos1.g 16000 0 1 && ./aos2.g 16000 0 1  &&  ./soa1.g 16000 0 1  && ./soa2.g 16000 0 1  && ./soa3.g 16000 0 1  && ./soa4.g 16000 0 1   && python3 error.py >> result_precision_one.txt && cat result_precision_one.txt	
	


run_init_512:
	./aos1.g 16000 0 0 && ./aos2.g 16000 0 0  &&  ./soa1.g 16000 0 0  && ./soa2.g 16000 0 0  && ./soa3.g 16000 0 0  && ./soa4.g 16000 0 0  && ./soa5.g 16000 0 0  && python3 error.py >> result_precision_init.txt && cat result_precision_init.txt

run_init_2:
	./aos1.g 16000 0 0 && ./aos2.g 16000 0 0  &&  ./soa1.g 16000 0 0  && ./soa2.g 16000 0 0  && ./soa3.g 16000 0 0  && ./soa4.g 16000 0 0  && python3 error.py >> result_precision_init.txt && cat result_precision_init.txt


