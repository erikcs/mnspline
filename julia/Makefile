# Makefile for OSX  (brew gcc) / Windows (with MinGW)
CFLAGS=-fopenmp -O3

ifeq ($(OS),Windows_NT)
	TARGET = mnspline.dll
	CC=gcc
endif
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	TARGET=mnspline.so
	CC=gcc-5
	CFLAGS += -fPIC
endif
# Linux...

mnspline:
	$(CC) ../src/mnspline.c -o $(TARGET) -shared $(CFLAGS)

