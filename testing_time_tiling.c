#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
//#include <cilk/cilk.h>
#include "xmmintrin.h"
#include "pmmintrin.h"
#include "space_orders_1.h"
#include "space_orders_1.h"
#include <unistd.h>

//#define min(a, b) (((a) < (b)) ? (a) : (b))

int main(int argc, char **argv)
{
  //int nrows = atoi(argv[1]);  //Number of rows
  //int ncols = atoi(argv[2]);  //Number of columns
  //int timesteps = atoi(argv[3]); //Number of timsteps
  int tile_size = atoi(argv[1]);   // Size of tile (space and time)
  int num_threads = atoi(argv[2]); // Size of tile (space and time)

  // Print a segment of the calculations and define starting and ending indexes
  int print_results = 0;
  int si = 1;
  int ei = 7;
  int sc = 1;
  int ec = 7;

  //Define sources and receivers numbers
  int nsrc = 1;   // Number of sources
  int nrecs = 10; // Number of receivers
  int omp_opt = 0; // Option for openMP . (NOT USED) TODO

  // Time measurement
  struct timeval t1, t2;
  double elapsedTime1, elapsedTime2;
  double sum_elapsedTime1, sum_elapsedTime2;

  omp_set_num_threads(num_threads);
  /* Flush denormal numbers to zero in hardware */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  //printf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n",nrows, ncols, timesteps, timestamps);

  //Define a 2d matrix to store the results
  double **RESULTS;
  malloc2d(&RESULTS, 40, 10);

  int validation_iters = 5; //Number of iterations
  int ri = 0;
  int rj = 0;
  int rs = atoi(argv[3]);
  int re = rs + 3;
  int ts = atoi(argv[4]);
  int te = ts + 3;

  for (int rows_pow = rs; rows_pow < re; rows_pow++)
  {
    for (int time_pow = ts; time_pow < te; time_pow++)
    {
      int nrows = pow(2, rows_pow);
      int ncols = pow(2, rows_pow);
      int timesteps = pow(2, time_pow);

      int timestamps = 2;  // To be removed after implementing buffering

      if (pow(2, 2 * rows_pow) * 2 < pow(2, 31))
      {
        

        RESULTS[ri][0] = nrows;
        RESULTS[ri][1] = timesteps;
        RESULTS[ri][2] = tile_size;

        printf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n", nrows, ncols, timesteps, timestamps);

        //printf("Allocated\n");

        sum_elapsedTime1 = 0;
        sum_elapsedTime2 = 0;
        double ***u;
        u = createMatrix(timestamps, nrows, ncols);
        double ***u2;
        u2 = createMatrix(timestamps, nrows, ncols);

        for (int validation_index = 0; validation_index < validation_iters; validation_index++)
        {
          // Set initial values

          initialize3(timestamps, nrows, ncols, u);
          initialize3(timestamps, nrows, ncols, u2);
          //printf("allocated");

          //printf("\n Starting Jacobi...");
          gettimeofday(&t1, NULL);
          u = so2_jacobi_3d(timesteps, nrows, ncols, u, omp_opt = 0, tile_size);
          gettimeofday(&t2, NULL);
          //printf("... Finished \n");
          elapsedTime1 = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000;
          //printf("Jacobi OpenMP, Time taken by program is : %3.3f\n",elapsedTime1);

          //sleep(2);
          //printf("Starting Jacobi...");
          gettimeofday(&t1, NULL);
          //u2 = so2_jacobi_3d(timesteps, nrows, ncols, u2, omp_opt = 0, tile_size);
          u2 = so2_tiled_skewed_jacobi_3d(timesteps, nrows, ncols, u2, omp_opt = 0, tile_size);
          gettimeofday(&t2, NULL);
          //printf("... Finished \n");
          elapsedTime2 = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000;
          //printf("Jacobi skewed, Time taken by program is : %3.3f\n",elapsedTime2);

          printf("Run : %d Speedup from skewing is : %3.3f\n", validation_index, elapsedTime1 / elapsedTime2);
          sum_elapsedTime1 += elapsedTime1;
          sum_elapsedTime2 += elapsedTime2;

          int validate_flag = 1;
          if (validate_flag)
          {
            for (int i = 0; i < nrows; i++)
            {
              for (int j = 0; j < ncols; j++)
              {
                if ((u[1][i][j] - u2[1][i][j]) > 0.00001)
                {
                  printf(" Failed %d, %d %f \n", i, j, (u[1][i][j] - u2[1][i][j]));
                }
              }
            }
          }
          if (print_results)
          {
            printf("\n ------------------------\n");
            for (int i = si; i < ei; i++)
            {
              printf("\n");
              for (int j = sc; j < ec; j++)
              {
                printf(" %3.3f", u[1][i][j]);
              }
            }
          }

          if (print_results)
          {
            printf("\n ------------------------\n");
            for (int i = si; i < ei; i++)
            {
              printf("\n");
              for (int j = sc; j < ec; j++)
              {
                printf(" %3.3f", u2[1][i][j]);
              }
            }
          }
        }

        RESULTS[ri][3] = sum_elapsedTime1 / validation_iters;
        RESULTS[ri][4] = sum_elapsedTime2 / validation_iters;
        RESULTS[ri][5] = sum_elapsedTime1 / sum_elapsedTime2;
        ri++;
        printf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n", nrows, ncols, timesteps, timestamps);
        printf(" Average speedup from skewing is : %3.3f\n", sum_elapsedTime1 / sum_elapsedTime2);

        free(u);
        free(u2);

        printf("\n");

      } //rows_pow
    } //if
  }   //time_pow

  if (1)
  {
    printf("\n ------------------------\n");
    for (int i = 0; i < ((re-rs)*(te-ts)); i++)
    {
      printf("\n");
      for (int j = 0; j < 6; j++)
      {
        printf(" %3.3f", RESULTS[i][j]);
      }
    }
  }

  char str[100] = "Results";
  // Uncomment to edit the name of the file.
  // printf("\n Enter the filename :"); gets(str);

  create_results_csv(str, RESULTS, 6, 20);
  return 0;
}
