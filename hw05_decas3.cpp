#include <stdio.h>  
#include <stdlib.h>  
#include <mpi.h>  
#include <queue>
  
#define iter_time 65536
#define MASTER_RANK 0  
#define TASK_TAG 142
#define RESULT_TAG 857
  
struct TaskData {  
    double currentX[4];
    double currentY[4];
    int len;
    double t;
};  
  
int main(int argc, char **argv) {  
    MPI_Init(&argc, &argv);  
  
    int size, rank;  
    MPI_Comm_size(MPI_COMM_WORLD, &size);  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  

    double step = 1.0 / iter_time;    
    double pointX[4], pointY[4];
    for (int i = 0; i < 4; i++) {  
        pointX[i] = i;  
        pointY[i] = i*i+1;  
    }
  
    if (rank == MASTER_RANK) {  
        // 主进程代码  
        std::queue<TaskData> q; // 模拟任务队列  

        for (double t = 0; t <= 1; t += step) {
            TaskData td;
            for (int i = 0; i < 4; i++) {
                td.currentX[i] = pointX[i];
                td.currentY[i] = pointY[i];
                td.len = 4;
                td.t = t;
            }
            q.push(td);
        }

        struct timespec sts,ets;
        timespec_get(&sts, TIME_UTC);
        
        while (!q.empty()){
            for (int worker = 1; worker < size; worker++) {
                // 发送任务给子进程  
                TaskData task_data;
                task_data = q.front(); q.pop();
                MPI_Send(&task_data, sizeof(TaskData), MPI_BYTE, worker, TASK_TAG, MPI_COMM_WORLD);  
            
                // 接收子进程的结果
                TaskData result;
                MPI_Recv(&result, sizeof(TaskData), MPI_BYTE, worker, RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (result.len > 1) {
                    q.push(result);
                }
                // printf("Master received from Rank %d: X = %f, Y = %f, len = %d\n", worker, result.currentX[0], result.currentY[0], result.len);
            } 
        }
        for (int worker = 1; worker < size; worker++) {
            TaskData end_data;
            end_data.len = -1;
            MPI_Send(&end_data, sizeof(TaskData), MPI_BYTE, worker, TASK_TAG, MPI_COMM_WORLD);
            // printf("send end to %d\n", worker);
        }
        timespec_get(&ets, TIME_UTC);
        time_t dsec=ets.tv_sec - sts.tv_sec;
        long dnsec=ets.tv_nsec - sts.tv_nsec;
        if (dnsec < 0){
            dsec--;
            dnsec += 1000000000ll;
        }
        printf("runtime: %ld.%09lds\n",dsec,dnsec);
    } else {  
        // 子进程代码  
        while (1) {
            TaskData task;  
            MPI_Recv(&task, sizeof(TaskData), MPI_BYTE, MASTER_RANK, TASK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
            
            if (task.len == -1) {  
                // printf("end!\n");
                break; // 退出循环  
            } 

            for (int j = 0; j < task.len-1; j++) {
                task.currentX[j] = (1-task.t) * task.currentX[j] + task.t * task.currentX[j+1];  
                task.currentY[j] = (1-task.t) * task.currentY[j] + task.t * task.currentY[j+1];  
            }
            task.len--;
            
            // 将结果发送回主进程  
            MPI_Send(&task, sizeof(TaskData), MPI_BYTE, MASTER_RANK, RESULT_TAG, MPI_COMM_WORLD);  
        }
    }  
  
    MPI_Finalize();  
    return 0;  
}