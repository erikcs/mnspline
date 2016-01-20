CC=clang-omp
CLFLAGS=-fopenmp -O3 -ffast-math

mnspline:
	$(CC) ../src/mnspline.c -o mnspline.so -shared -fPIC $(CLFLAGS)
