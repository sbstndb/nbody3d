all: add_mul add pearson  obj

add_mul: add_mul.c
	gcc -g -march=native -mtune=native -Ofast -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize -mavx $< -o $@ -lm -lSDL2
	#clang -g -march=native -mtune=native -O3 -fno-math-errno -funroll-loops -finline-functions -ftree-vectorize -mavx $< -o $@ -lm -lSDL2

add: add.c
	gcc -g -march=native -mtune=native -Ofast -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize -mavx $< -o $@ -lm -lSDL2
	#clang -g -march=native -mtune=native -O3 -fno-math-errno -funroll-loops -finline-functions -ftree-vectorize -mavx $< -o $@ -lm -lSDL2



pearson: pearson.c
	gcc -g -march=native -mtune=native -O3 -funroll-loops -finline-functions -fpeel-loops -ftree-vectorize -ftree-loop-vectorize -mavx $< -o $@ -lm -lSDL2
	#clang -g -march=native -mtune=native -O3 -fno-math-errno -funroll-loops -finline-functions -ftree-vectorize -mavx $< -o $@ -lm -lSDL2




obj:
	objdump -d add_mul >> add_mul.s
	objdump -d add >> add.s
	objdump -d pearson >> pearson.s

clean:
	rm -Rf *~ add_mul add pearson
