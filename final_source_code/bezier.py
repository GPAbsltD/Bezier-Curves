import taichi as ti
import time
import timeit

ti.init(arch=ti.cpu)

iter_time = 16
# 定义点的坐标为二维向量场
points = ti.Vector.field(2, dtype=ti.f32, shape=4)  # 每个元素存储一个二维坐标(x, y)
# 结果存储向量场，每个元素对应一个时间步的坐标(x, y)
result = ti.Vector.field(2, dtype=ti.f32, shape=iter_time+1)


@ti.kernel
def bernstein_parallel():
    step = 1.0 / iter_time
    points_mat = ti.Matrix.zero(ti.f32, 4, 2)

    # 初始化点坐标
    for i in range(points.shape[0]):
        points_mat[i, :] = [i, i * i + 1]

    # 并行循环计算并存储结果
    for t in range(int(iter_time + 1)):
        t_float = t * step

        # current = (1 - t_float) ** 3 * points[0] + t_float * (1 - t_float) ** 2 * points[1] + \
        #           (t_float ** 2) * (1 - t_float) * points[2] + (t_float ** 3) * points[3]

        coeff = ti.Matrix([[ti.pow((1 - t_float), 3),
                            3 * t_float * ti.pow((1 - t_float), 2),
                            3 * ti.pow(t_float, 2) * (1 - t_float),
                            ti.pow(t_float, 3)]])
        current = coeff @ points_mat
        result[t] = current[0, :]


if __name__ == '__main__':
    start = time.time()
    bernstein_parallel()
    end = time.time()
    # print("runtime: ", end-start, sep='')

    exec_time = timeit.timeit(bernstein_parallel, number=1000)
    # print(exec_time)

    exec_time = timeit.timeit(bernstein_parallel, number=10000)
    # print(exec_time)

    # 打印结果
    for t in range(iter_time+1):
        print('t=', t/iter_time, '  X=', result[t][0], '  Y=', result[t][1], sep='')
