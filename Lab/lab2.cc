#include <mpi.h>
#include <stdio.h>
#include <cstdio>
#include <math.h>

int main(int argc, char** argv){
	MPI_Init(&argc,&argv);
	int rank,size;
	double temp_sum= 0;
	double final_sum= 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	int k = atoi(argv[1]);

	for(int i =rank;i<k;i+=size){
	 	temp_sum += 4* sqrt((1-(i/(double)k)*(i/(double)k)))/(double)k;
	}
	

	//MPI_Reduce(&sendbuf,&recvbuf,count,datatype,op,comm)
	MPI_Allreduce(&temp_sum,&final_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		if(rank==0){
			printf("The total sum of area is %f \n",final_sum);
		}
	
	
	MPI_Finalize();
	return 0;
}



