#define _GNU_SOURCE
#include <sched.h>
#include <thread>
#include <iostream>
#include <cmath>
#include <string>

void cal(int rank, int cpu, double slice, double *global){
    double local = 0;
    for(double i = rank; i <slice; i +=cpu){
        local += sqrt(1.0 - (i/slice)*(i/slice)) / slice ;
    }
    *global = local;

    /*std::cout << rank << '  ';
    std::cout << local << std::endl;*/
}

int main(int argc, char* argv[]){
    double slice =  atoi(argv[1]);

    cpu_set_t cpuset;
    sched_getaffinity(0, sizeof(cpuset), &cpuset);
    int cpu = CPU_COUNT(&cpuset);

    double global[cpu], ans=0;

    std::thread sread[cpu];
    for (int i = 0; i < cpu; i++) {
        sread[i] = std::thread(cal, i, cpu, slice, &global[i]);
    }
    for(int i=0; i<cpu; i++){
        sread[i].join();
    }
    for(int i=0; i<cpu; i++){
        ans += global[i];
    }

    printf("%.6lf\n",ans*4);

    return 0;
    }
