tile_size = 16;  %Size of tile (space and time)

%Print a segment of the calculations and define starting and ending indexes
print_results = 0;
si = 1;
ei = 7;
sc = 1;
ec = 7;

%Define sources and receivers numbers
nsrc = 1;   %Number of sources
nrecs = 10; % Number of receivers
omp_opt = 0; % Option for openMP . (NOT USED) TODO

%Time measurement
%struct timeval t1, t2;
%double elapsedTime1, elapsedTime2;
%double sum_elapsedTime1, sum_elapsedTime2;


fprintf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n",nrows, ncols, timesteps, timestamps);

%Define a 2d matrix to store the results

Results = zeros(40,10);

validation_iters = 5; %Number of iterations
ri = 0;
rj = 0;
rs = 6;
re = rs + 3;
ts = 6;
te = ts + 3;


for rows_pow = rs:re
  for time_pow = ts:te
      nrows = 2^rows_pow;
      ncols = 2^rows_pow;
      timesteps = 2^time_pow;
      timestamps = 2;  %To be removed after implementing buffering

      if 2^(2 * rows_pow) * 2 < (2^31)
        Results(ri,0) = nrows;
        Results(ri,1) = timesteps;
        Results(ri,2) = tile_size;

        fprintf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n", nrows, ncols, timesteps, timestamps);
        %printf("Allocated\n");
        sum_elapsedTime1 = 0;
        sum_elapsedTime2 = 0;
        
        u(1:timestamps,1:nrows,1:ncols) = 0;
        u2(1:timestamps,1:nrows,1:ncols) = 0;
        
        
        src_coords(1:nsrc,1:2) = 0;
        src(1:timesteps,1:nsrc) = 0;

        for isd= (1:2:timesteps)
          src(isd,0) = isd^2 - isd;
        end

        src_coords(0,0) = floor(nrows/2);
        src_coords(0,1) = floor(ncols/2);

        for validation_index = (0:1:validation_iters)
          %Set initial values

          initialize3(timestamps, nrows, ncols, u);
          initialize3(timestamps, nrows, ncols, u2);
          %printf("allocated");

          %printf("\n Starting Jacobi...");
          tic;
          %u = jacobi_3d_src_rcv(timesteps, nrows, ncols, u, src_coords, src, omp_opt=0, tile_size);
          toc;
          %printf("... Finished \n");
          %elapsedTime1 = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000;
          %printf("Jacobi OpenMP, Time taken by program is : %3.3f\n",elapsedTime1);

          %sleep(2);
          %printf("Starting Jacobi...");
          tic;
          %u2 = jacobi_3d_validate(timesteps, nrows, ncols, u2, src_coords, src, omp_opt=0, tile_size);
          toc;
          %printf("... Finished \n");
          %elapsedTime2 = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec) / 1000000;
          %printf("Jacobi skewed, Time taken by program is : %3.3f\n",elapsedTime2);

          fprintf("Run : %d Speedup from skewing is : %3.3f\n", validation_index, elapsedTime1 / elapsedTime2);
          %sum_elapsedTime1 += elapsedTime1;
          %sum_elapsedTime2 += elapsedTime2;

          validate_flag = 1;
          
          if validate_flag == 1
            for i = (0:1:nrows)
              for j = (0:1:ncols)
                if ((u(1,i,j) - u2(1,i,j) > 0.00001))
                  fprintf(" Failed %d, %d %f \n", i, j, (u(1,i,j) - u2(1,i,j)));
                
                end
              end
            end
          end
          if print_results==1
          
            printf("\n ------------------------\n");
            for i = (si:1:ei)
              printf("\n");
              for j = (sc:1:ec)
                 printf(" %3.3f", u(1,i,j));
              end
            end
          end

          if print_results==1
          
            fprintf("\n ------------------------\n");
            for i = (si:1:ei)
              printf("\n");
              for j = (sc:1:ec)
                  printf(" %3.3f", u2(1,i,j));
              end
            end
          end
        end

        Results(ri,3) = sum_elapsedTime1 / validation_iters;
        Results(ri,4) = sum_elapsedTime2 / validation_iters;
        Results(ri,5) = sum_elapsedTime1 / sum_elapsedTime2;
        ri = ri +1;
        fprintf("Problem setup is \nnrows: %d\nncols: %d\nTimesteps: %d\nTimestamps: %d\n", nrows, ncols, timesteps, timestamps);
        fprintf(" Average speedup from skewing is : %3.3f\n", sum_elapsedTime1 / sum_elapsedTime2);

        %free(u);
        %free(u2);

        printf("\n");

     end %rows_pow
  end %if
end   %time_pow

  if 1
    printf("\n ------------------------\n");
    for i = (0:1:((re-rs)*(te-ts)))
      printf("\n");
      for j = (0:1:6)
        printf(" %3.3f", Results(i,j));
      end 
    end
  end

  char str[100] = "Results";
  create_results_csv(str, Results, 6, 20);
 