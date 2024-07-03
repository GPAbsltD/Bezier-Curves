#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
using namespace std;

#define iter_time 65536
int main(int argc, char **argv) {  
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    double step = 1.0 / iter_time;
    double t_start = 1.0 * rank / size;
    double t_end = (rank == size - 1) ? 1.0 : 1.0 * (rank+1) / size - step;

    // printf("Rank %d: t_start = %f, t_end = %f\n", rank, t_start, t_end);  
 
    double pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {  
        pointX[i] = i;  
        pointY[i] = i*i+1;  
    }

    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);
  
    for (double t = t_start; t < t_end; t += step) {  
        double currentX, currentY;  
        currentX = pow((1-t), 3)*pointX[0] + t*pow((1-t), 2)*pointX[1] +
                    pow(t, 2)*(1-t)*pointX[2] + pow(t, 3)*pointX[3];
        currentY = pow((1-t), 3)*pointY[0] + t*pow((1-t), 2)*pointY[1] +
                    pow(t, 2)*(1-t)*pointY[2] + pow(t, 3)*pointY[3];
        // draw(currentX, currentY), 忽略在实验中  
        // printf("Rank %d: t = %f, X = %f, Y = %f\n", rank, t, currentX, currentY);  
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