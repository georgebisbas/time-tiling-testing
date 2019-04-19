# time-tiling-testing
A repository to test loop nest optimizations and particularly time-tiling

Run `make`
and ./a.out `tile_size` `number of threads` `size as power of 2` `timesteps as power of 2` (e.g. ./a.out 32 8 8 5 (for `tile_size=32` , 8 threads, square grid of 2^arg2 per dimension, timesteps = 2^arg3 ))

at the end of the execution you can name the csv file containing the results.(Default name is Results.csv)

For thread pinning:

use
For `icc` compiler
```KMP_AFFINITY="proclist=[0-7]omp" ./a.out arg0 8 arg2 arg3```

```KMP_AFFINITY="proclist=[0-3]omp" ./a.out arg0 4 arg2 arg3```

For `gcc` compiler
```GOMP_CPU_AFFINITY="0-7" ./a.out arg0 8 arg2 arg3```

```GOMP_CPU_AFFINITY="0-3" ./a.out arg0 4 arg2 arg3``` 

Currently experiments are executed for user-defined parameters and scale up to 16*grid and 4*timesteps... to be updated soon...
