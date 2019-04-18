#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#include <stdio.h>
#include <string.h>

// Function to malloc a 2d array
double ***so2_jacobi_3d(int timesteps, int nrows, int ncols, double ***grid, int omp_opt, int tile_size)
{
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1 = 1;
  int t0 = 0;

  int x_m = 2;
  int x_M = nrows - 2;
  int y_m = 2;
  int y_M = ncols - 2;

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
        grid[t1][xi][yi] = (grid[t0][xi][yi] + grid[t0][xi - 1][yi] + grid[t0][xi + 1][yi] + grid[t0][xi][yi - 1] + grid[t0][xi][yi + 1]+
                                               grid[t0][xi - 2][yi] + grid[t0][xi + 2][yi] + grid[t0][xi][yi - 2] + grid[t0][xi][yi + 2]) / 5; }
    }
  }
  return grid;
}

double ***so2_tiled_skewed_buffered(int timesteps, int nrows, int ncols, double ***grid, int omp_opt, int tile_size)
{
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1 = 1;
  int t0 = 0;

  int x_m = 2;
  int x_M = nrows - 2;
  int y_m = 2;
  int y_M = ncols - 2;

  int t_blk_size = tile_size;
  int x_blk_size = tile_size;

  // x blocking
  for (int xblk = x_m; xblk < x_M + timesteps - 2; xblk += x_blk_size)
  {
    // t interchange
    for (int titer = 2; titer < timesteps - 2; titer++)
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
          grid[t1][xi - titer][yi - titer] = (grid[t0][xi - titer][yi - titer] + grid[t0][(xi - titer) - 1][yi - titer] + grid[t0][(xi - titer) + 1][yi - titer] + grid[t0][xi - titer][yi - titer - 1] + grid[t0][xi - titer][yi - titer + 1]+ grid[t0][xi - 1][yi] + grid[t0][xi + 1][yi] + grid[t0][xi][yi -1] + grid[t0][xi][yi + 1]) / 5;
        }
      }
    }
  }
  return grid;
}


double ***so2_tiled_skewed_jacobi_3d(int timesteps, int nrows, int ncols, double ***grid, int omp_opt, int tile_size)
{
  int t = 0;
  int xi = 0;
  int yi = 0;
  int t1 = 1;
  int t0 = 0;

  int x_m = 2;
  int x_M = nrows - 2;
  int y_m = 2;
  int y_M = ncols - 2;

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
          //printf(" Speedup from skeweing is : u[%d][%d][%d]\n", titer,xi,yi);

          grid[t1][xi - titer][yi - titer] = (grid[t0][xi - titer][yi - titer] + grid[t0][(xi - titer) - 1][yi - titer] + grid[t0][(xi - titer) + 1][yi - titer] + grid[t0][xi - titer][yi - titer - 1] + grid[t0][xi - titer][yi - titer + 1] + grid[t0][(xi - titer) - 2][yi - titer] + grid[t0][(xi - titer) + 2][yi - titer] + grid[t0][xi - titer][yi - titer - 2] + grid[t0][xi - titer][yi - titer + 2])/ 5;}
      }
    }
  }
  return grid;
}
