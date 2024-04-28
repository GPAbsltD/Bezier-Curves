#include <iostream>
using namespace std;

void deCas() {
    double currentX[4], currentY[4];
    const int iter_time = 1024;
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

#include <math.h>
void Bern(){
    const int iter_time = 1024;
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

#include <emmintrin.h>

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

void Bern_para(){
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

    const int iter_time = 1024;
    float step = 1.0 / iter_time;
    
    for (float t = 0; t <= 1.0; t = t + step) {
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

/*
// x86 运算时间计算 精度不够
#include <sys/time.h>
#include <cstdlib>

int main() {
    clock_t start, finish;
    start = clock();
    Bern_para();
    finish = clock();

    cout << "serial runtime: " << endl;
    cout << (finish-start)/float(CLOCKS_PER_SEC) << " s " << endl;
    return 0;
}
*/

#include <time.h>

int main() {
    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);
    Bern_para();
    timespec_get(&ets, TIME_UTC);

    time_t dsec=ets.tv_sec - sts.tv_sec;
    long dnsec=ets.tv_nsec - sts.tv_nsec;
    if (dnsec < 0){
        dsec--;
        dnsec += 1000000000ll;
    }
    printf("decas runtime: %ld.%09lds\n",dsec,dnsec);
    return 0;
}
