#include <stdio.h>
#include <stdlib.h>
#include <thread>

using namespace std;
const int INF = ((1 << 30) - 1);
const int V = 50010;
void input(char* inFileName, int rank, int num_cpu);
void output(char* outFileName/*, int rank, int num_cpu*/);

void block_FW(int B, int num_cpu);
int ceil(int a, int b);
void cal(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height, int rank, int num_cpu);

int n, m;
static int Dist[V][V];

int main(int argc, char* argv[]) {

    int num_cpu;
    cpu_set_t cpuset; 
    sched_getaffinity(0, sizeof(cpuset), &cpuset);
    num_cpu = CPU_COUNT(&cpuset); 
    thread th[num_cpu];
    for (int i = 0; i < num_cpu; i++)
    {
        th[i] = thread(input, argv[1], i, num_cpu);
    }
    for(int i = 0; i < num_cpu; i++){
        th[i].join();
    }
    //input(argv[1]);
    int B = 512;
    block_FW(B, num_cpu);

/*    thread th0[num_cpu];
    for (int i = 0; i < num_cpu; i++)
    {
        th0[i] = thread(output, argv[2], i, num_cpu);
    }
    for(int i = 0; i < num_cpu; i++){
        th0[i].join();
    }*/
    output(argv[2]);
    return 0;
}

inline void input(char* infile, int rank, int num_cpu) {
    FILE* file = fopen(infile, "rb");
    fread(&n, sizeof(int), 1, file);
    fread(&m, sizeof(int), 1, file);
#pragma unroll
    for (int i = rank; i < n; i+= num_cpu) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                Dist[i][j] = 0;
            } else {
                Dist[i][j] = INF;
            }
        }
    }

    int pair[3];
    for (int i = 0; i < m; ++i) {
        fread(pair, sizeof(int), 3, file);
        Dist[pair[0]][pair[1]] = pair[2];
    }
    fclose(file);
}

/*inline void output(char* outFileName, int rank, int num_cpu) {
    FILE* outfile = fopen(outFileName, "w");
#pragma unroll
    for (int i = rank; i < n; i+= num_cpu) {
        for (int j = 0; j < n; ++j) {
            if (Dist[i][j] >= INF) Dist[i][j] = INF;
        }
        fwrite(Dist[i], sizeof(int), n, outfile);
    }
    fclose(outfile);
}*/
void output(char* outFileName) {
    FILE* outfile = fopen(outFileName, "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (Dist[i][j] >= INF) Dist[i][j] = INF;
        }
        fwrite(Dist[i], sizeof(int), n, outfile);
    }
    fclose(outfile);
}

int ceil(int a, int b) { return (a + b - 1) / b; }

inline void block_FW(int B, int num_cpu) {
    int round = ceil(n, B);
    /*find the # of cpu */

#pragma unroll
    for (int r = 0; r < round; ++r) {
        
        /* Phase 1*/
        thread th1[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th1[i] = thread(cal, B, r, r, r, 1, 1, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th1[i].join();
        }

        /* Phase 2*/
        thread th2[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th2[i] = thread(cal, B, r, r, 0, r, 1, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th2[i].join();
        }

        thread th3[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th3[i] = thread(cal, B, r, r, r + 1, round - r - 1, 1, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th3[i].join();
        }

        thread th4[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th4[i] = thread(cal, B, r, 0, r, 1, r, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th4[i].join();
        }

        thread th5[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th5[i] = thread(cal, B, r, r + 1, r, 1, round - r - 1, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th5[i].join();
        }
        // cal(B, r, r, r + 1, round - r - 1, 1);
        // cal(B, r, 0, r, 1, r);
        // cal(B, r, r + 1, r, 1, round - r - 1);

        /* Phase 3*/
        thread th6[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th6[i] = thread(cal, B, r, 0, 0, r, r, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th6[i].join();
        }

        thread th7[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th7[i] = thread(cal, B, r, 0, r + 1, round - r - 1, r, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th7[i].join();
        }

        thread th8[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th8[i] = thread(cal, B, r, r + 1, 0, r, round - r - 1, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th8[i].join();
        }

        thread th9[num_cpu];
        for(int i = 0; i < num_cpu; ++i){
            th9[i] = thread(cal, B, r, r + 1, r + 1, round - r - 1, round - r - 1, i, num_cpu);
        }
        for(int i = 0; i < num_cpu; ++i){
            th9[i].join();
        }        
        // cal(B, r, 0, 0, r, r);
        // cal(B, r, 0, r + 1, round - r - 1, r);
        // cal(B, r, r + 1, 0, r, round - r - 1);
        // cal(B, r, r + 1, r + 1, round - r - 1, round - r - 1);
    }
}

inline void cal(
    int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height, int rank, int num_cpu) {
    int block_end_x = block_start_x + block_height;
    int block_end_y = block_start_y + block_width;
#pragma unroll
    for (int b_i = block_start_x+ rank; b_i < block_end_x; b_i += num_cpu) {
        for (int b_j = block_start_y ; b_j < block_end_y; ++b_j) {
            // To calculate B*B elements in the block (b_i, b_j)
            // For each block, it need to compute B times
            for (int k = Round * B; k < (Round + 1) * B && k < n; ++k) {
                // To calculate original index of elements in the block (b_i, b_j)
                // For instance, original index of (0,0) in block (1,2) is (2,5) for V=6,B=2
                int block_internal_start_x = b_i * B;
                int block_internal_end_x = (b_i + 1) * B;
                int block_internal_start_y = b_j * B;
                int block_internal_end_y = (b_j + 1) * B;

                if (block_internal_end_x > n) block_internal_end_x = n;
                if (block_internal_end_y > n) block_internal_end_y = n;

                for (int i = block_internal_start_x; i < block_internal_end_x; ++i) {
                    for (int j = block_internal_start_y; j < block_internal_end_y; ++j) {
                        if (Dist[i][k] + Dist[k][j] < Dist[i][j]) {
                            Dist[i][j] = Dist[i][k] + Dist[k][j];
                        }
                    }
                }
            }
        }
    }
}
