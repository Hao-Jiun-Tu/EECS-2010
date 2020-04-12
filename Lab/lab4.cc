#include <inttypes.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
int64_t neg_for_composite(int64_t n) {
    for (int64_t i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            return -n;
        }
    }
    return n;
}
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int64_t n = atoll(argv[1]);
    int64_t partial_sum = 0;
    int64_t total_sum = 0;
    #pragma omp parallel for reduction(+:partial_sum)
        for(int64_t i = rank + 2;i <= n; i += size){
            partial_sum += neg_for_composite(i);
        }
        MPI_Reduce(&partial_sum, &total_sum, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(rank == 0){
        printf("%" PRId64 "\n", total_sum);
    }
    
    MPI_Finalize();
}
