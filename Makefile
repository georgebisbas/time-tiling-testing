all: BUFFER_time_tiling

CC = gcc #icc
DEPS = tile_auxiliary.h


BUFFER_time_tiling: BUFFER_time_tiling.c
	$(CC) -O3 -g -ffast-math -w -fopenmp BUFFER_time_tiling.c -lm
clean:
	rm -f BUFFER_time_tiling *.o
