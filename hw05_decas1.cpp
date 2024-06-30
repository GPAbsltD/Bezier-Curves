#include <stdio.h>  
#include <stdlib.h>  
#include <mpi.h>  
#include <time.h>
  
#define iter_time 65536
  
int main(int argc, char **argv) {  
    MPI_Init(&argc, &argv);  
  
    int size, rank;  
    MPI_Comm_size(MPI_COMM_WORLD, &size);  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
  
    double step = 1.0 / iter_time;
    double t_start = 1.0 * rank / size;
    double t_end = (rank == size - 1) ? 1.0 : 1.0 * (rank+1) / size - step;

    double pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {
        pointX[i] = i;
        pointY[i] = i*i+1;
    }

    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);

    double currentX[4], currentY[4];
    for (double t = t_start; t < t_end; t += step) {
        for (int i = 0; i < 4; i++) {
            currentX[i] = pointX[i];
            currentY[i] = pointY[i];
        }

        for (int i = 3; i > 0; i--) {
            for (int j = 0; j < i; j++) {
                currentX[j] = (1-t) * currentX[j] + t * currentX[j+1];  
                currentY[j] = (1-t) * currentY[j] + t * currentY[j+1];  
            }  
        }
        // printf("Rank %d: t = %f, X = %f, Y = %f\n", rank, t, currentX[0], currentY[0]);
    }  
    
    timespec_get(&ets, TIME_UTC);
    time_t dsec=ets.tv_sec - sts.tv_sec;
    long dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("runtime: %ld.%09lds\n",dsec,dnsec);

    MPI_Finalize();  
    return 0;  
}