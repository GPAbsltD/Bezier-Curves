#include <stdio.h>  
#include <stdlib.h>  
#include <pthread.h>  
#include <math.h>
#include <omp.h>
  
#define iter_time 524288
#define MAX_THREADS 8  
  
// 线程参数结构  
typedef struct {   
    double pointX[4];  
    double pointY[4];  
    int thread_id;
    double step;
} ThreadArgs;  
  
void* computeBezierPoint(void* args_void) {  
    ThreadArgs* args = (ThreadArgs*)args_void;  
    int t_id = args->thread_id;
    double pointX[4] = {args->pointX[0], args->pointX[1], args->pointX[2], args->pointX[3]};  
    double pointY[4] = {args->pointY[0], args->pointY[1], args->pointY[2], args->pointY[3]};  
    float step = args->step;
  
    for (float t = float(t_id) / MAX_THREADS; t < float(t_id+1) / MAX_THREADS; t += step) {
        double currentX[4], currentY[4];
        for (int i = 0; i < 4; i++) {
            currentX[i] = pointX[i];
            currentY[i] = pointY[i];
        }

        for (int i = 3; i > 0; i--) {
            for (int j = 0; j < i; j++){
                currentX[j] = (1-t) * currentX[j] + t * currentX[j+1];
                currentY[j] = (1-t) * currentY[j] + t * currentY[j+1];
            }
        }

        // draw (currentX[0], currentY[0])
	}      
    return NULL;  
}  
  
void deCas_pth() {    
    // test cases  
    double pointX[4], pointY[4];  
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

void deCas_omp(){
    double step = 1.0 / iter_time;

    // test cases
    double pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {
        pointX[i] = i;
        pointY[i] = i*i+1;
    }
    
    omp_set_num_threads(MAX_THREADS);
    #pragma omp parallel for
	for (int i = 0; i < iter_time; i++) {
        double t = double(i)/iter_time;
        double currentX[4], currentY[4];
        for (int i = 0; i < 4; i++) {
            currentX[i] = pointX[i];
            currentY[i] = pointY[i];
        }

        for (int i = 3; i > 0; i--) {
            for (int j = 0; j < i; j++){
                currentX[j] = (1-t) * currentX[j] + t * currentX[j+1];
                currentY[j] = (1-t) * currentY[j] + t * currentY[j+1];
            }
        }

        // draw (currentX[0], currentY[0])
	}
}

void deCas() {
    double currentX[4], currentY[4];
    double step = 1.0 / iter_time;

    // test cases
    double pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {
        pointX[i] = i;
        pointY[i] = i*i+1;
    }
    
	for (float t = 0; t <= 1.0; t = t + step) {
        for (int i = 0; i < 4; i++) {
            currentX[i] = pointX[i];
            currentY[i] = pointY[i];
        }

        for (int i = 3; i > 0; i--) {
            for (int j = 0; j < i; j++){
                currentX[j] = (1-t) * currentX[j] + t * currentX[j+1];
                currentY[j] = (1-t) * currentY[j] + t * currentY[j+1];
            }
        }
        // draw(currentX[0], currentY[0]), ignoring in the experiment
	}
}
  
int main() {
    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);
    deCas_pth();  
    timespec_get(&ets, TIME_UTC);

    time_t dsec=ets.tv_sec - sts.tv_sec;
    long dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("Pthread runtime: %ld.%09lds\n",dsec,dnsec);

    timespec_get(&sts, TIME_UTC);
    deCas_omp();  
    timespec_get(&ets, TIME_UTC);

    dsec=ets.tv_sec - sts.tv_sec;
    dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("OpenMP runtime: %ld.%09lds\n",dsec,dnsec);

    timespec_get(&sts, TIME_UTC);
    deCas();  
    timespec_get(&ets, TIME_UTC);

    dsec=ets.tv_sec - sts.tv_sec;
    dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("Serial runtime: %ld.%09lds\n",dsec,dnsec);
    return 0;  
}