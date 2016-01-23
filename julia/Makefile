# Makefile for OSX / Windows (with MinGW)

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	TARGET = mnspline.so
	CC=clang-omp
else
	TARGET = mnspline.dll
	CC=gcc
endif

CLFLAGS=-fopenmp -O3

mnspline:
	$(CC) ../src/mnspline.c -o $(TARGET) -shared -fPIC $(CLFLAGS)

