#include <stdio.h>  
#include <stdlib.h>  
#include <pthread.h>  
#include <math.h>
#include <omp.h>
  
#define iter_time 65536
#define MAX_THREADS 8  
  
// 线程参数结构  
typedef struct {   
    double pointX[4];  
    double pointY[4];  
    int thread_id;
    double step;
    // 如果需要输出或其他处理，可以在这里添加额外的字段  
} ThreadArgs;  
  
// 线程函数  
void* computeBezierPoint(void* args_void) {  
    ThreadArgs* args = (ThreadArgs*)args_void;  
    int t_id = args->thread_id;
    double pointX[4] = {args->pointX[0], args->pointX[1], args->pointX[2], args->pointX[3]};  
    double pointY[4] = {args->pointY[0], args->pointY[1], args->pointY[2], args->pointY[3]};  
    float step = args->step;
  
    for (float t = float(t_id) / MAX_THREADS; t < float(t_id+1) / MAX_THREADS; t += step) {
        float currentX =    pow((1-t), 3)*pointX[0]     + t*pow((1-t), 2)*pointX[1] + 
                            pow(t, 2)*(1-t)*pointX[2]   + pow(t, 3)*pointX[3];
        float currentY =    pow((1-t), 3)*pointY[0]     + t*pow((1-t), 2)*pointY[1] + 
                            pow(t, 2)*(1-t)*pointY[2]   + pow(t, 3)*pointY[3];
	}

    // draw (currentX, currentY)
  
    return NULL;  
}  
  
void Bern_pth() {    
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


void Bern_omp(){
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
        float currentX =    pow((1-t), 3)*pointX[0]     + t*pow((1-t), 2)*pointX[1] + 
                            pow(t, 2)*(1-t)*pointX[2]   + pow(t, 3)*pointX[3];
        float currentY =    pow((1-t), 3)*pointY[0]     + t*pow((1-t), 2)*pointY[1] + 
                            pow(t, 2)*(1-t)*pointY[2]   + pow(t, 3)*pointY[3];
        // draw(currentX, currentY), ignoring in the experiment
	}
}


void Bern(){
    double step = 1.0 / iter_time;

    // test cases
    double pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {
        pointX[i] = i;
        pointY[i] = i*i+1;
    }
    
	for (float t = 0; t <= 1.0; t = t + step) {
        float currentX =    pow((1-t), 3)*pointX[0]     + t*pow((1-t), 2)*pointX[1] + 
                            pow(t, 2)*(1-t)*pointX[2]   + pow(t, 3)*pointX[3];
        float currentY =    pow((1-t), 3)*pointY[0]     + t*pow((1-t), 2)*pointY[1] + 
                            pow(t, 2)*(1-t)*pointY[2]   + pow(t, 3)*pointY[3];
        // draw(currentX, currentY), ignoring in the experiment
	}
}
  
int main() {
    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);
    Bern_pth();  
    timespec_get(&ets, TIME_UTC);

    time_t dsec=ets.tv_sec - sts.tv_sec;
    long dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("Pthread runtime: %ld.%09lds\n",dsec,dnsec);

    timespec_get(&sts, TIME_UTC);
    Bern_omp();  
    timespec_get(&ets, TIME_UTC);

    dsec=ets.tv_sec - sts.tv_sec;
    dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("OpenMP runtime: %ld.%09lds\n",dsec,dnsec);

    timespec_get(&sts, TIME_UTC);
    Bern();  
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