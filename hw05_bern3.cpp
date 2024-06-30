#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <emmintrin.h>
using namespace std;

#define iter_time 100
#define MAX_THREADS 8

float sse_dot_product(__m128 a, __m128 b) {  
    // 对应元素相乘  
    __m128 mul_result = _mm_mul_ps(a, b);  
    
    float dot_product[4];
    _mm_store_ps(dot_product, mul_result);

    float summ = 0;
    for (int i = 0; i < 4; i++) {
        summ += dot_product[i];
    }

    return summ;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int start_i = iter_time * rank / size;
    int end_i = iter_time * (rank+1) / size - 1;
  
    float array[4];
    __m128 row0 = _mm_set_ps(-1.0f, 3.0f, -3.0f, 1.0f);
    __m128 row1 = _mm_set_ps(0.0f, 3.0f, -6.0f, 3.0f);
    __m128 row2 = _mm_set_ps(0.0f, 0.0f, 3.0f, -3.0f);
    __m128 row3 = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
    
    // 加载控制点X和Y到SSE寄存器
    float pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {
        pointX[i] = i;
        pointY[i] = i*i+1;
    }
    __m128 controlX = _mm_load_ps(pointX);
    __m128 controlY = _mm_load_ps(pointY);

    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);
    
    omp_set_num_threads(MAX_THREADS);
    #pragma omp parallel for
	for (int i = start_i; i <= end_i; i++) {
        double t = double(i) / iter_time;
        __m128 T = _mm_set_ps(1, t, t*t, t*t*t);

        // 先算后半部分的乘法
        float rowT[4];
        rowT[0] = sse_dot_product(row0, T);
        rowT[1] = sse_dot_product(row1, T);
        rowT[2] = sse_dot_product(row2, T);
        rowT[3] = sse_dot_product(row3, T);
        
        __m128 rowT_vec = _mm_load_ps(rowT);

        // 提取结果到单独的变量中 
        float resultX = sse_dot_product(controlX, rowT_vec);
        float resultY = sse_dot_product(controlY, rowT_vec);

        // printf("Rank %d: t = %f, X = %f, Y = %f\n", rank, t, resultX, resultY);
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