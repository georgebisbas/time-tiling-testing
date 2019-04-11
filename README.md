# time-tiling-testing
A repository to test loop nest optimizations and particularly time-tiling

Run `make`
and ./a.out `tile_size` `number of threads` (e.g. ./a.out 32 8 (for `tile_size=32` and 8 threads))

at the end of the execution you will be asked for a name for the csv file containing the results.

Warning : do not change the problem parameters until time buffering is implemented. Soon to be.

For thread pinning:

use
For `icc` compiler
```make && KMP_AFFINITY="proclist=[0-7]omp" ./a.out 32 8```
```make && KMP_AFFINITY="proclist=[0-3]omp" ./a.out 32 4```

For `gcc` compiler
```make && GOMP_CPU_AFFINITY="0-7" ./a.out 32 8```
```make && GOMP_CPU_AFFINITY="0-3" ./a.out 32 4``` 
