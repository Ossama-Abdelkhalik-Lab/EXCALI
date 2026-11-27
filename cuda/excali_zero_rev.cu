// CUDA kernel for EXCALI's zero-revolution path (prototype, compiled but never run:
// no GPU was available). One thread per case; inputs in Fortran (3,n) layout.
#include "excali_zero_rev.cuh"
#include <cuda_runtime.h>

extern "C" __global__ void excali_zero_rev_kernel(int n, const double *__restrict__ r1, const double *__restrict__ r2,
                                                  const double *__restrict__ tof, const int *__restrict__ dir,
                                                  double *__restrict__ v1, double *__restrict__ v2,
                                                  int *__restrict__ status) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double a[3] = {r1[3 * i], r1[3 * i + 1], r1[3 * i + 2]};
    double b[3] = {r2[3 * i], r2[3 * i + 1], r2[3 * i + 2]};
    double w1[3], w2[3];
    int st, half;
    exc::excali_zero_rev(a, b, tof[i], dir[i], w1, w2, st, half);
    for (int k = 0; k < 3; ++k) { v1[3 * i + k] = w1[k]; v2[3 * i + k] = w2[k]; }
    status[i] = st;
}

// Host wrapper: copies n cases to the device, runs the kernel, copies back.
// Returns the CUDA error code (0 on success).
extern "C" int excali_zero_rev_gpu(int n, const double *r1, const double *r2, const double *tof, const int *dir,
                                   double *v1, double *v2, int *status) {
    double *d_r1, *d_r2, *d_t, *d_v1, *d_v2; int *d_dir, *d_st;
    size_t s3 = 3 * sizeof(double) * (size_t)n, s1 = sizeof(double) * (size_t)n, si = sizeof(int) * (size_t)n;
    cudaError_t e;
    if ((e = cudaMalloc(&d_r1, s3)) || (e = cudaMalloc(&d_r2, s3)) || (e = cudaMalloc(&d_t, s1)) ||
        (e = cudaMalloc(&d_dir, si)) || (e = cudaMalloc(&d_v1, s3)) || (e = cudaMalloc(&d_v2, s3)) ||
        (e = cudaMalloc(&d_st, si))) return (int)e;
    cudaMemcpy(d_r1, r1, s3, cudaMemcpyHostToDevice); cudaMemcpy(d_r2, r2, s3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_t, tof, s1, cudaMemcpyHostToDevice); cudaMemcpy(d_dir, dir, si, cudaMemcpyHostToDevice);
    int bs = 128;
    excali_zero_rev_kernel<<<(n + bs - 1) / bs, bs>>>(n, d_r1, d_r2, d_t, d_dir, d_v1, d_v2, d_st);
    e = cudaGetLastError();
    if (!e) e = cudaMemcpy(v1, d_v1, s3, cudaMemcpyDeviceToHost);
    if (!e) e = cudaMemcpy(v2, d_v2, s3, cudaMemcpyDeviceToHost);
    if (!e) e = cudaMemcpy(status, d_st, si, cudaMemcpyDeviceToHost);
    cudaFree(d_r1); cudaFree(d_r2); cudaFree(d_t); cudaFree(d_dir); cudaFree(d_v1); cudaFree(d_v2); cudaFree(d_st);
    return (int)e;
}
