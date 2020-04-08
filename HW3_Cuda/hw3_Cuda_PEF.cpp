#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <iostream>

#define DEV_NO 0
#define ROUND_MAX 4

const int INF = ((1 << 30) - 1);
inline void input(char *inFileName);
inline void output(char *outFileName);
void block_FW(int B);
int ceil(int a, int b);
__global__ void Phase1_Kernel(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height, int n, int* Dist_gpu, int pitch_int);
__global__ void Phase2_Kernel(int B, int Round, int n, int* Dist_gpu, int pitch_int);
__global__ void Phase3_Kernel(int B, int Round, int n, int* Dist_gpu, int pitch_int);

int n, m;
int* Dist;
int* Dist_gpu;
cudaDeviceProp prop;
size_t pitch;

int main(int argc, char* argv[])
{
    cudaGetDeviceProperties(&prop, DEV_NO);
    input(argv[1]);
    int B = 64;
    output(argv[2]);
    return 0;
}

inline void input(char* infile) {
    FILE* file = fopen(infile, "rb");
    fread(&n, sizeof(int), 1, file);
    fread(&m, sizeof(int), 1, file);
    cudaMallocHost( &Dist, sizeof(int) * n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) Dist[i * n + j] = 0;
            else Dist[i * n + j] = INF;
        }
    }
    int pair[3];
    for (int i = 0; i < m; ++i) {
        fread(pair, sizeof(int), 3, file);
        Dist[pair[0] * n + pair[1]] = pair[2];
    }
    fclose(file);
}

void output(char *outFileName) {
    FILE *outfile = fopen(outFileName, "w");
    fwrite(Dist, sizeof(int), n * n, outfile);
    fclose(outfile);
}

int ceil(int a, int b) {
    return (a + b -1)/b;
}

void block_FW(int B)
{
    int thread_num = (B * B > prop.maxThreadsPerBlock)? prop.maxThreadsPerBlock: B * B;
    int round = (n + B - 1) / B;
    cudaMallocPitch((void**) &Dist_gpu, &pitch, sizeof(int) * n, n + 64);
    int pitch_int = pitch / sizeof(int);
    cudaMemcpy2D(Dist_gpu, pitch, Dist, sizeof(int) * n, sizeof(int) * n, n, cudaMemcpyHostToDevice);

    dim3 grid2(round - 1, 2);
    dim3 grid3(round - 1, round - 1);
    dim3 block(B, thread_num / B);
    for (int r = 0; r < round; ++r) {
        /* Phase 1*/
        Phase1_Kernel<<<1, block, B * B * sizeof(int) >>>(B, r, r, r, 1, 1, n, Dist_gpu, pitch_int);
        /* Phase 2*/
        Phase2_Kernel<<<grid2, block, 2 * B * B * sizeof(int)>>>(B, r, n, Dist_gpu, pitch_int);
        // cal_kernel<<<r*1, thread_num>>>(B, r, r, 0, r, 1, n, Dist_gpu);
        // cal_kernel<<<(round - r - 1), thread_num>>>(B, r, r, r + 1, round - r -1, 1, n, Dist_gpu);
        // cal_kernel<<<r, thread_num>>>(B, r, 0, r, 1, r, n, Dist_gpu);
        // cal_kernel<<< 1*(round-r-1), thread_num>>>(B, r, r + 1, r, 1, round - r -1, n, Dist_gpu);
        /* Phase 3*/
        Phase3_Kernel<<<grid3, block>>>(B, r, n, Dist_gpu, pitch_int);
        // cal_kernel<<<r * r, thread_num>>>(B, r, 0, 0, r, r, n, Dist_gpu);
        // cal_kernel<<<(round - r - 1) * r, thread_num>>>(B, r, 0, r + 1, round - r - 1, r, n, Dist_gpu);
        // cal_kernel<<<r * (round - r - 1), thread_num>>>(B, r, r + 1, 0, r, round - r - 1, n, Dist_gpu);
        // cal_kernel<<< (round - r - 1) * (round - r - 1), thread_num>>>(B, r, r + 1, r + 1, round - r - 1, round - r - 1, n, Dist_gpu);
    }
    cudaMemcpy2D(Dist, sizeof(int) * n, Dist_gpu, pitch, sizeof(int) * n, n, cudaMemcpyDeviceToHost);
}

__global__ void Phase1_Kernel(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height, int n, int* Dist_gpu, int pitch_int) {
    int b_i = block_start_x ;
    int b_j = block_start_y ;

    extern __shared__ int shared_mem[];
    int global_i[ROUND_MAX],
        global_j[ROUND_MAX],
        inner_i[ROUND_MAX],
        inner_j[ROUND_MAX];

    #pragma unroll
    for(int r = 0; r < 4; r++){
        inner_i[r] = threadIdx.y + 16 * r;
        inner_j[r] = threadIdx.x;

        global_i[r] = b_i * B + inner_i[r];
        global_j[r] = b_j * B + inner_j[r];
        if (!((global_i[r] >= n) | (global_j[r] >= n)))
            shared_mem[inner_i[r] * B + inner_j[r]] = Dist_gpu[global_i[r] * pitch_int + global_j[r]];
    }
    for (int k = 0; k <  B && (k + Round * B) < n; ++k) {
        __syncthreads();
        #pragma unroll
        for(int r = 0; r < 4; r++){
            if ((global_i[r] >= n) | (global_j[r] >= n)) continue;
            if (shared_mem[inner_i[r] * B + inner_j[r]] > shared_mem[inner_i[r] * B + k] + shared_mem[k * B + inner_j[r]])
                shared_mem[inner_i[r] * B + inner_j[r]] = shared_mem[inner_i[r] * B + k] + shared_mem[k * B + inner_j[r]];
        }
    }
    #pragma unroll
    for(int r = 0; r < 4; r++){
        if (!((global_i[r] >= n) | (global_j[r] >= n)))
            Dist_gpu[global_i[r] * pitch_int + global_j[r]] = shared_mem[inner_i[r] * B + inner_j[r]];
    }
}

__global__ void Phase2_Kernel(int B, int Round, int n, int* Dist_gpu, int pitch_int) {
    int b_i, b_j;
    extern __shared__ int shared_mem[];
    if(blockIdx.y == 0){
        b_i = Round;
        b_j = blockIdx.x + (blockIdx.x >= Round);
    }
    else{
        b_i = blockIdx.x + (blockIdx.x >= Round);
        b_j = Round;
    }
    int global_i[ROUND_MAX],
        global_j[ROUND_MAX],
        inner_i[ROUND_MAX],
        inner_j[ROUND_MAX];

    #pragma unroll
    for(int r = 0; r < 4; r++){
        inner_i[r] = threadIdx.y + 16 * r;
        inner_j[r] = threadIdx.x;
        global_i[r] = b_i * B + inner_i[r];
        global_j[r] = b_j * B + inner_j[r];
        int global_pivot_i = Round * B + inner_i[r];
        int global_pivot_j = Round * B + inner_j[r];
        if (!((global_i[r] >= n) | (global_j[r] >= n)))
            shared_mem[inner_i[r] * B + inner_j[r]] = Dist_gpu[global_i[r] * pitch_int + global_j[r]];
        if (!((global_pivot_i >= n) | (global_pivot_j >= n)))
            shared_mem[inner_i[r] * B + inner_j[r] + B * B] = Dist_gpu[global_pivot_i * pitch_int + global_pivot_j];
    }

    for (int k = 0; k < B && (k + Round * B) < n; ++k) {
        __syncthreads();
        #pragma unroll
        for(int r = 0; r < 4; r++){
            if ((global_i[r] >= n) | (global_j[r] >= n)) continue ;
            if (shared_mem[inner_i[r] * B + inner_j[r]] > shared_mem[inner_i[r] * B + k + !blockIdx.y * B * B] + shared_mem[k * B + inner_j[r] + blockIdx.y * B * B])
                shared_mem[inner_i[r] * B+inner_j[r]] = shared_mem[inner_i[r] * B + k + !blockIdx.y * B * B] + shared_mem[k * B+inner_j[r] + blockIdx.y * B * B];
        }
    }
    #pragma unroll
    for(int r = 0; r < 4; r++){
        if (!((global_i[r] >= n) | (global_j[r] >= n)))
            Dist_gpu[global_i[r] * pitch_int + global_j[r]] = shared_mem[inner_i[r] * B + inner_j[r]];
    }
}

__global__ void Phase3_Kernel(int B, int Round, int n, int* Dist_gpu, int pitch_int) {
    int b_i = blockIdx.y + (blockIdx.y >= Round);
    int b_j = blockIdx.x + (blockIdx.x >= Round);

    __shared__ int shared_mem[8192];

    int global_i[ROUND_MAX],
        global_j[ROUND_MAX],
        inner_i[ROUND_MAX],
        inner_j[ROUND_MAX];

    int my_dist[ROUND_MAX];

    #pragma unroll
    for(int r = 0; r < 4; r++){
        inner_i[r] = threadIdx.y + 16 * r;
        inner_j[r] = threadIdx.x;
        global_i[r] = b_i * B + inner_i[r];
        global_j[r] = b_j * B + inner_j[r];
        int row_pivot_i = global_i[r];
        int row_pivot_j = Round * B + inner_j[r];
        int col_pivot_i = Round * B + inner_i[r];
        int col_pivot_j = global_j[r];

        my_dist[r] = Dist_gpu[global_i[r] * pitch_int + global_j[r]];
        shared_mem[inner_i[r] * B + inner_j[r]] = Dist_gpu[row_pivot_i * pitch_int + row_pivot_j];
        shared_mem[inner_i[r] * B + inner_j[r] + B * B] = Dist_gpu[col_pivot_i * pitch_int + col_pivot_j];
    }
    __syncthreads();
    for (int k = 0; k <  B && (k + Round * B) < n; ++k) {
        #pragma unroll
        for(int r = 0; r < 4; r++){
            int tmp = shared_mem[inner_i[r] * B + k] + shared_mem[k * B + inner_j[r] + B * B];
            if (my_dist[r] > tmp)  my_dist[r] = tmp;
        }
    }
    #pragma unroll
    for(int r = 0; r < 4; r++)
        Dist_gpu[global_i[r] * pitch_int + global_j[r]] = my_dist[r];
}

