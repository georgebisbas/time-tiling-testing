#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#include <stdio.h>
#include <string.h>

//Function to store results in a csv file
void create_results_csv(char *filename, double **a, int n, int m)
{
  printf("\n Creating %s.csv file", filename);
  FILE *fp;
  int i, j;
  filename = strcat(filename, ".csv");
  fp = fopen(filename, "w+");
  fprintf(fp, "NROWS, NCOLS, Timesteps, Tilesize, el1, el2, speedup");
  for (i = 0; i < m; i++)
  {
    fprintf(fp, "\n%d", i + 1);
    for (j = 0; j < n; j++)
      fprintf(fp, ", %f ", a[i][j]);
  }
  fclose(fp);
  printf("\n %sfile created", filename);
}

// Function to allocate memory for u[t][x][y]
double ***createMatrix(int timestamps, int nrows, int ncols)
{
  //  printf("\nAllocating 3d...");
  double ***matrix;

  matrix = (double ***)malloc(timestamps * sizeof(double **));
  for (int nt = 0; nt < timestamps; nt++)
  {
    matrix[nt] = (double **)malloc(nrows * sizeof(double *));
    for (int row = 0; row < nrows; row++)
    {
      matrix[nt][row] = (double *)malloc(ncols * sizeof(double));
    }
  }
  //printf("...DONE\n");
  return matrix;
}

// Function to allocate memory for u[t][x][y]
void initialize3(int timestamps, int nrows, int ncols, double ***A_init)
{
  // Initialize the matrix
  //printf("\nInit 3d...");
  int ti = 0;
  int xi = 0;
  int yi = 0;
  for (ti = 0; ti < timestamps; ti += 1)
  {
    for (xi = 0; xi < nrows; xi += 1)
    {
      for (yi = 0; yi < ncols; yi += 1)
      {
        A_init[ti][xi][yi] = 2 * xi - 1.4 * yi + 0.5 * fabs(xi - yi);
      }
    }
  }
  //printf("...DONE\n");
}

// Function to malloc a 2d array
double malloc2d(double ***C, int nrows, int ncols)
{
  int i;
  *C = malloc(sizeof(double *) * nrows);

  if (*C == NULL)
  {
    printf("ERROR: out of memory\n");
    return 1;
  }

  for (i = 0; i < nrows; i++)
  {
    (*C)[i] = malloc(sizeof(double) * ncols);
    if ((*C)[i] == NULL)
    {
      printf("ERROR: out of memory\n");
      return 1;
    }
  }
  //printf("Allocated!\n");
  return 0;
}

// Function to malloc a 2d array
double ***jacobi_3d(int timesteps, int nrows, int ncols, double ***grid, int omp_opt, int tile_size)
{
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1 = 1;
  int t0 = 0;
  int x_m = 1;
  int x_M = nrows - 1;
  int y_m = 1;
  int y_M = ncols - 1;

  for (int titer = 0; titer < timesteps - 1; titer++)
  {
    t1 = titer % 2;
    t0 = !t1;

#pragma omp parallel for
    for (int xi = x_m; xi < x_M; xi++)
    {
#pragma omp simd
      for (int yi = y_m; yi < y_M; yi++)
      {
        grid[t1][xi][yi] = (grid[t0][xi][yi] + grid[t0][xi - 1][yi] + grid[t0][xi + 1][yi] + grid[t0][xi][yi - 1] + grid[t0][xi][yi + 1]) / 5;
      }
    }
  }
  return grid;
}

double ***tiled_skewed_jacobi_3d(int timesteps, int nrows, int ncols, double ***grid, int omp_opt, int tile_size)
{
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1 = 1;
  int t0 = 0;

  int x_m = 1;
  int x_M = nrows - 1;
  int y_m = 1;
  int y_M = ncols - 1;

  int t_blk_size = tile_size;
  int x_blk_size = tile_size;

  // x blocking
  for (int xblk = x_m; xblk < x_M + timesteps - 1; xblk += x_blk_size)
  {
    // t interchange
    for (int titer = 0; titer < timesteps - 1; titer++)
    {

      t1 = titer % 2;
      t0 = !t1;
// parallelize for x blocks
#pragma omp parallel for
      for (int xi = max(x_m + titer, xblk); xi < min((x_M + titer), (xblk + x_blk_size)); xi++)
      {
        // simd
#pragma omp simd
        for (int yi = y_m + titer; yi < y_M + titer; yi++)
        {
          grid[t1][xi - titer][yi - titer] = (grid[t0][xi - titer][yi - titer] + grid[t0][(xi - titer) - 1][yi - titer] + grid[t0][(xi - titer) + 1][yi - titer] + grid[t0][xi - titer][yi - titer - 1] + grid[t0][xi - titer][yi - titer + 1]) / 5;
        }
      }
    }
  }
  return grid;
}

double ***tiled_skewed_buffered(int timesteps, int nrows, int ncols, double ***grid, int omp_opt, int tile_size)
{
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1 = 1;
  int t0 = 0;

  int x_m = 1;
  int x_M = nrows - 1;
  int y_m = 1;
  int y_M = ncols - 1;

  int t_blk_size = tile_size;
  int x_blk_size = tile_size;

  // x blocking
  for (int xblk = x_m; xblk < x_M + timesteps - 1; xblk += x_blk_size)
  {
    // t interchange
    for (int titer = 0; titer < timesteps - 1; titer++)
    {
      t1 = titer % 2;
      t0 = !t1;
// parallelize for x blocks
#pragma omp parallel for
      for (int xi = max(x_m + titer, xblk); xi < min((x_M + titer), (xblk + x_blk_size)); xi++)
      {
        // simd
#pragma omp simd
        for (int yi = y_m + titer; yi < y_M + titer; yi++)
        {
          grid[t1][xi - titer][yi - titer] = (grid[t0][xi - titer][yi - titer] + grid[t0][(xi - titer) - 1][yi - titer] + grid[t0][(xi - titer) + 1][yi - titer] + grid[t0][xi - titer][yi - titer - 1] + grid[t0][xi - titer][yi - titer + 1]) / 5;
        }
      }
    }
  }
  return grid;
}