
all: nbody.g nbody1.g nbody2.g nbody3.g nbody3_clang.g #nbody.i

nbody.g: nbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

nbody1.g: nbody1.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize -mavx2 -Ofast -fopt-info-all=nbody1.gcc.optrpt $< -o $@ -lm -fopenmp

nbody2.g: nbody2.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize -mavx2 -Ofast -fopt-info-all=nbody2.gcc.optrpt $< -o $@ -lm -fopenmp


nbody3.g: nbody3.c
	gcc -march=native -mtune=native -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize -mavx2 -Ofast -fopt-info-all=nbody3.gcc.optrpt $< -o $@ -lm -fopenmp

nbody3_clang.g: nbody3.c
	clang -march=native -mtune=native -funroll-loops -finline-functions -ftree-vectorize -mavx2 -Ofast  $< -o $@ -lm -fopenmp

nbody.i: nbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

clean:
	rm -Rf *~ nbody.g nbody.i nbody1.g nbody2.g nbody3.g  nbody3_clang.g *.optrpt
	
run:
	./nbody.g && ./nbody1.g && ./nbody2.g && ./nbody3.g

