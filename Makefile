all: time_tiling

CC = gcc #icc
DEPS = tile_auxiliary.h


time_tiling: time_tiling.c
	$(CC) -O3 -g -ffast-math -w -fopenmp time_tiling.c -lm
clean:
	rm -f time_tiling *.o
