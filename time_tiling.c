#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cilk/cilk.h>
#include "xmmintrin.h"
#include "pmmintrin.h"
#include "tiled_buffer.h"

//#define min(a, b) (((a) < (b)) ? (a) : (b))





int main(int argc, char **argv) {
    //int nrows = atoi(argv[1]);  //Number of rows
    //int ncols = atoi(argv[2]);  //Number of columns
    //int timesteps = atoi(argv[3]); //Number of timsteps
    int tile_size = atoi(argv[1]); // Size of tile (space and time)
    int print_results = 0;    // Print a segment of the calculations
    int si = 10; int ei = 15; int sc = 10; int ec = 15; // Bounds of the segment

    int nsrc = 1;  // Number of sources
    int nrecs = 10; // Number of receivers
    //int timestamps = timesteps;     // To be removed after implementing buffering
    //int timestamps2 = timesteps;    // To be removed after implementing buffering
    int omp_opt=0;                  // Option for openMP . NOT USED

    // Time measurement
    struct timeval t1, t2;
    double elapsedTime1, elapsedTime2;
		double sum_elapsedTime1, sum_elapsedTime2;

    //printf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n",nrows, ncols, timesteps, timestamps);

    //double ** A;
    //malloc2d(&A, nrows, ncols);
    //double ** B; // B is like A(t-1)
    //malloc2d(&B, nrows, ncols);
    //double ** src_coords; // B is like A(t-1)
    //malloc2d(&src_coords, nsrc, 2);
    //double ** rec_coords; // B is like A(t-1)
    //malloc2d(&rec_coords, nsrc, 2);
    //double ** src; // B is like A(t-1)
    //malloc2d(&src, nsrc, 2);

		double ** RESULTS;
    malloc2d(&RESULTS, 40, 10);

		int validation_iters = 10;
		int ri=-1; int rj=0;
		for(int rows_pow = 9 ; rows_pow < 13; rows_pow++){
			//for(int cols_pow = 5 ; cols_pow< 11; cols_pow++){
//ri++;
			for(int time_pow = 2 ; time_pow < 6; time_pow++){
ri++;


			int nrows = pow( 2 ,  rows_pow);
			int ncols = pow( 2 ,  rows_pow);
			int timesteps = pow( 2 ,  time_pow);
			int timestamps = timesteps;     // To be removed after implementing buffering
			int timestamps2 = timesteps;    // To be removed after implementing buffering

			RESULTS[ri][0] = nrows;
			RESULTS[ri][1] = timesteps;
			RESULTS[ri][2] = tile_size;

		printf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n", nrows, ncols, timesteps, timestamps);

    double *** u;
    u = createMatrix(nrows, ncols, timestamps);
    double *** u2;
    u2 = createMatrix(nrows, ncols, timestamps2);

    /* Flush denormal numbers to zero in hardware */
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

		sum_elapsedTime1 = 0 ;
		sum_elapsedTime2 = 0;

for(int validation_index=0; validation_index<validation_iters;validation_index++){
    initialize3(nrows, ncols,timestamps, u);
    initialize3(nrows, ncols, timestamps2, u2);

    //printf("\n Starting Jacobi...");
		gettimeofday(&t1, NULL);
    u = buffer_jacobi_3d(timesteps, nrows, ncols, u, omp_opt=0, tile_size);
    gettimeofday(&t2, NULL);
    //printf("... Finished \n");
    elapsedTime1 = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    //printf("Jacobi OpenMP, Time taken by program is : %3.3f\n",elapsedTime1);

    //printf("Starting Jacobi...");
		gettimeofday(&t1, NULL);
    u2 = tiled_skewed_buffer_jacobi_3d(timesteps, nrows, ncols, u2, omp_opt=0,tile_size);
    gettimeofday(&t2, NULL);
    //printf("... Finished \n");
    elapsedTime2 = (double)(t2.tv_sec-t1.tv_sec)+(double)(t2.tv_usec-t1.tv_usec)/1000000;
    //printf("Jacobi skewed, Time taken by program is : %3.3f\n",elapsedTime2);

    printf(" Speedup from skeweing is : %3.3f\n", elapsedTime1/elapsedTime2);
		sum_elapsedTime1 += elapsedTime1;
		sum_elapsedTime2 += elapsedTime2;
}

RESULTS[ri][3] = sum_elapsedTime1/validation_iters;
RESULTS[ri][4] = sum_elapsedTime2/validation_iters;
RESULTS[ri][5] = sum_elapsedTime1/sum_elapsedTime2;



printf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n",nrows, ncols, timesteps, timestamps);
printf(" Average speedup from skeweing is : %3.3f\n", sum_elapsedTime1/sum_elapsedTime2);
//RESULTS[5][5]= sum_elapsedTime1/sum_elapsedTime2;



int validate_flag=1;
if(validate_flag){
    //for (int k = 0; k < timestamps; k++) {
        for (int i = 0; i < nrows; i++) {
          for (int j = 0; j < ncols; j++) {
          if ((u[i][j][timesteps-1] - u2[i][j][timesteps-1])> 0.001) {
            printf(" Failed %d, %d %f \n",i, j,(u[i][j][timesteps-1] - u2[i][j][timesteps-1]));
          }
        }
      }
    //}
}


if(print_results){
  printf("\n ------------------------\n");
    for (int i = si; i < ei; i++) {
      printf("\n");
      for (int j = sc; j < ec; j++) {
      printf(" %3.3f", u[i][j][timesteps-1]);
    }
  }
}

if(print_results){
  printf("\n ------------------------\n");
    for (int i = si; i < ei; i++) {
      printf("\n");
      for (int j = sc; j < ec; j++) {
      printf(" %3.3f", u2[i][j][timesteps-1]);
    }
  }
}


		free(u);
    free(u2);
    printf("\n");

	}//rows_pow
//}//cols_pow

}//time_pow

if(1){
  printf("\n ------------------------\n");
    for (int i = 0; i < 20; i++) {
      printf("\n");
      for (int j = 0; j < 6; j++) {
      printf(" %3.3f", RESULTS[i][j]);
    }
  }
}
    return 0;
}
