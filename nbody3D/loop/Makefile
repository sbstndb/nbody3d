all : 0 l v lv vl

lv:
	gcc loop.c -g -march=native -mtune=native -funroll-loops $(avx512) -mavx2 -mavx -ftree-vectorize -fno-math-errno -Ofast -lm -o lv

vl:
	gcc loop.c -g -march=native -mtune=native $(avx512) -mavx2 -mavx -ftree-vectorize -funroll-loops -fno-math-errno  -Ofast -lm -o vl

v:
	gcc loop.c -g -march=native -mtune=native $(avx512) -mavx2 -mavx -ftree-vectorize -fno-unroll-loops -fno-math-errno  -Ofast -lm -o v

l:
	gcc loop.c -g -march=native -mtune=native $(avx512) -mavx2 -mavx -fno-tree-vectorize -funroll-loops -funroll-all-loops -fno-math-errno  -Ofast -lm -o l

0:
	gcc loop.c -g -march=native -mtune=native $(avx512) -mavx2 -mavx -fno-tree-vectorize -fno-unroll-loops -fno-math-errno  -Ofast -lm -o 0

	
clean:
	rm l v lv vl 0
	
run:
	make all
	./0 ./l ./v ./lv ./vl
	
