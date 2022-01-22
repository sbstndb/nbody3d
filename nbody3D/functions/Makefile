all : vectorized non_vectorized

vectorized:
	gcc benchmark_functions.c -march=native -mtune=native $(avx512) -mavx2 -mavx -ftree-vectorize -funroll-loops -O3 -lm -o v3
	gcc benchmark_functions.c -march=native -mtune=native $(avx512) -mavx2 -mavx -ftree-vectorize -funroll-loops -Ofast -lm -o vf

non_vectorized:
	gcc benchmark_functions.c -march=native -mtune=native -fno-tree-vectorize -funroll-loops -O3 -lm -o nv3
	gcc benchmark_functions.c -march=native -mtune=native -fno-tree-vectorize -funroll-loops -Ofast -lm -o nvf
	
clean:
	rm v3 vf nv3 nvf
	
run:
	make all
	./nf3 ./nvf ./v3 ./vf
	
