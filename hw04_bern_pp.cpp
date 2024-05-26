#include <stdio.h>  
#include <stdlib.h>  
#include <pthread.h>  
#include <math.h>
#include <omp.h>
#include <emmintrin.h>
#include <iostream>
  
#define iter_time 524288
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

void omp_pp(){
    // 加载矩阵的行到SSE寄存器  
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

    float step = 1.0 / iter_time;
    
    omp_set_num_threads(MAX_THREADS);
    #pragma omp parallel for
	for (int i = 0; i < iter_time; i++) {
        double t = double(i)/iter_time;
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
	}
}
  
// 线程参数结构  
typedef struct {   
    float pointX[4];  
    float pointY[4];  
    int thread_id;
    double step;
    // 如果需要输出或其他处理，可以在这里添加额外的字段  
} ThreadArgs;  
  
// 线程函数  
void* computeBezierPoint(void* args_void) {  
    ThreadArgs* args = (ThreadArgs*)args_void;  
    int t_id = args->thread_id;
    float pointX[4] = {args->pointX[0], args->pointX[1], args->pointX[2], args->pointX[3]};  
    float pointY[4] = {args->pointY[0], args->pointY[1], args->pointY[2], args->pointY[3]};  
    float step = args->step;

    __m128 row0 = _mm_set_ps(-1.0f, 3.0f, -3.0f, 1.0f);  
    __m128 row1 = _mm_set_ps(0.0f, 3.0f, -6.0f, 3.0f);  
    __m128 row2 = _mm_set_ps(0.0f, 0.0f, 3.0f, -3.0f);  
    __m128 row3 = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);  
    __m128 controlX = _mm_load_ps(pointX);
    __m128 controlY = _mm_load_ps(pointY);
  
    for (float t = float(t_id) / MAX_THREADS; t < float(t_id+1) / MAX_THREADS; t += step) {
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
	}

    // draw (currentX, currentY)
  
    return NULL;  
}  
  
void pth_pp() {    
    // test cases  
    float pointX[4], pointY[4];  
    for (int i = 0; i < 4; i++) {  
        pointX[i] = i;  
        pointY[i] = i*i+1;  
    }  
  
    // 线程数组和ID  
    pthread_t threads[MAX_THREADS];  
    ThreadArgs* args[MAX_THREADS];  
    float step = float(1) / iter_time;
  
    // 为每个t值创建线程  
    for (int i = 0; i < MAX_THREADS; i++) {  
        // 分配并初始化线程参数  
        args[i] = (ThreadArgs*)malloc(sizeof(ThreadArgs));  
        args[i]->thread_id = i;  
        for (int j = 0; j < 4; j++) {  
            args[i]->pointX[j] = pointX[j];  
            args[i]->pointY[j] = pointY[j];  
        }  
        args[i]->step = step;
  
        // 创建线程  
        pthread_create(&threads[i], NULL, computeBezierPoint, (void*)args[i]);  
    }  
  
    // 等待所有线程完成  
    for (int i = 0; i < MAX_THREADS; i++) {  
        pthread_join(threads[i], NULL);  
    }  
}  
  
int main() {
    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);
    pth_pp();  
    timespec_get(&ets, TIME_UTC);

    time_t dsec=ets.tv_sec - sts.tv_sec;
    long dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("Pthread runtime: %ld.%09lds\n",dsec,dnsec);

    timespec_get(&sts, TIME_UTC);
    omp_pp();  
    timespec_get(&ets, TIME_UTC);

    dsec=ets.tv_sec - sts.tv_sec;
    dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("OpenMP runtime: %ld.%09lds\n",dsec,dnsec);
    return 0;  
}