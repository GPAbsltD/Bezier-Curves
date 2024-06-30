#include <stdio.h>  
#include <stdlib.h>  
#include <mpi.h>
#include <omp.h>
  
#define iter_time 100
#define MAX_THREADS 8
  
int main(int argc, char **argv) {  
    MPI_Init(&argc, &argv);  
  
    int size, rank;  
    MPI_Comm_size(MPI_COMM_WORLD, &size);  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
  
    int start_i = iter_time * rank / size;
    int end_i = iter_time * (rank+1) / size - 1;

    double pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {
        pointX[i] = i;
        pointY[i] = i*i+1;
    }

    double currentX[4], currentY[4];
    omp_set_num_threads(MAX_THREADS);
    #pragma omp parallel for
    for (int i = start_i; i <= end_i; i++) {  
        double t = double(i) / iter_time;
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
        printf("Rank %d: t = %f, X = %f, Y = %f\n", rank, t, currentX[0], currentY[0]);
    }  
  
    MPI_Finalize();  
    return 0;  
}