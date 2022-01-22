
all: aos1.g aos2.g soa1.g soa2.g soa3.g soa4.g soa5.g

aos1.g:aos1.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp
	
	
aos2.g: aos2.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -Ofast $< -o $@ -lm -fopenmp
	

soa1.g:soa1.c
	gcc -march=native -mtune=native -fno-unroll-loops -fno-inline-functions -fno-peel-loops -fno-tree-vectorize -fno-tree-loop-vectorize -Ofast $< -o $@ -lm -fopenmp
	
soa2.g:soa2.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp
	
soa3.g:soa3.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp
	
soa4.g:soa4.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp			

soa5.g:soa5.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize $(avx) -mavx2 -fno-math-errno -Ofast $< -o $@ -lm -fopenmp


aos1.i: aos1.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

clean:
	rm -Rf *~ *.g *.i *.optrpt *.txt
	
run:
	./aos1.g && ./aos2.g &&  ./soa1.g && ./soa2.g && ./soa3.g && ./soa4.g && ./soa5.g

