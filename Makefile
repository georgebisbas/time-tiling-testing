all: testing_time_tiling

CC = gcc #icc
DEPS = tile_auxiliary.h


testing_time_tiling: testing_time_tiling.c
	$(CC) -O3 -g -ffast-math -w -fopenmp testing_time_tiling.c -lm
clean:
	rm -f testing_time_tiling *.o
