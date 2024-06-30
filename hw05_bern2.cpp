#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
using namespace std;

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
    
    omp_set_num_threads(MAX_THREADS);
    #pragma omp parallel for
    for (int i = start_i; i <= end_i; i++) {  
        double t = double(i) / iter_time;
        double currentX, currentY;  
        currentX = pow((1-t), 3)*pointX[0] + t*pow((1-t), 2)*pointX[1] +
                    pow(t, 2)*(1-t)*pointX[2] + pow(t, 3)*pointX[3];
        currentY = pow((1-t), 3)*pointY[0] + t*pow((1-t), 2)*pointY[1] +
                    pow(t, 2)*(1-t)*pointY[2] + pow(t, 3)*pointY[3];
        // draw(currentX, currentY), 忽略在实验中  
        printf("Rank %d: t = %f, X = %f, Y = %f\n", rank, t, currentX, currentY);  
    }
  
    MPI_Finalize();  
    return 0;  
}