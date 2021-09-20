CC=g++
CFLAGS=-O3 -fopenmp -march=native -std=c++17 -I include

ifeq ($(shell uname -s), Darwin)
	CC=clang++
	CFLAGS=-Xpreprocessor -fopenmp -lomp -O3 -march=native -std=c++11 -I include
endif

main:
	$(CC) main.cpp -o loopyBP $(CFLAGS)
