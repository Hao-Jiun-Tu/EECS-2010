#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <float.h>
#include <cstring>
#include <algorithm>
#include <boost/sort/pdqsort/pdqsort.hpp>
using namespace std;

int main(int argc, char** argv){

        int rank, size;
        int n = atoi(argv[1]);
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_File f;
        MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
        int procs_data_num;
        float *cd;
        if(size > n){
                if(rank == 0){
                    cd = new float [n];
		    procs_data_num = n;
                    MPI_File_read_at(f, 0, cd, n, MPI_FLOAT, MPI_STATUS_IGNORE);
		    boost::sort::pdqsort(cd, cd + procs_data_num);
                }
                else
                    MPI_File_read_at(f, 0, cd, 0, MPI_FLOAT, MPI_STATUS_IGNORE);
        }
        else if(n % size != 0){
                if(rank < n % size){
                    procs_data_num = n / size + 1;
                    cd = new float [procs_data_num];
                    MPI_File_read_at(f, sizeof(float) * (procs_data_num) * rank, cd, procs_data_num, MPI_FLOAT, MPI_STATUS_IGNORE);
                }
                else if(rank < size){
                    procs_data_num = n / size;
                    cd = new float [procs_data_num];
                    MPI_File_read_at(f, sizeof(float)*procs_data_num*rank + sizeof(float)*(n%size), cd, procs_data_num, MPI_FLOAT, MPI_STATUS_IGNORE);
                }
        }
        else{
            procs_data_num = n / size;
	    cd = new float [procs_data_num];
            MPI_File_read_at(f, sizeof(float) * (procs_data_num) * rank, cd, (procs_data_num), MPI_FLOAT, MPI_STATUS_IGNORE);
        }

	if(size <= n) boost::sort::pdqsort(cd, cd + procs_data_num);

	float *reg;
        float *reg_2;

	if(size > n){
	}
       	else if(n % size == 0){
                reg = new float [procs_data_num];
                reg_2 = new float [procs_data_num];
            for (int i = 0;i < size / 2 + 1 ; i++) {
                //Even Phase
                if(rank % 2 == 0 && rank != size - 1){
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                        //Even Rank Sent Message To Odd Rank
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //Even Rank Recv Message From Odd Rank
                        int i, j, k;
                        for(i = j = k = 0; k < procs_data_num; k++){
				if(cd[i] > reg[j]) reg_2[k] = reg[j++];
				else reg_2[k] = cd[i++];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 1){
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank - 1, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);
                        int i, j, k;
                        for(i = j = k = procs_data_num - 1; k >= 0; k--){
				if(cd[i] >= reg[j]) reg_2[k] = cd[i--];
				else   reg_2[k] = reg[j--];
			}
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
		else{}

                //Odd Phase
                if(rank % 2 == 1 && rank != size - 1){
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                        //Odd Rank Sent Message To Even Rank
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //Odd Rank Recv Message From Even Rank
                        int i, j, k;
                        for(i = j = k = 0; k < procs_data_num; k++){
                            if(cd[i] > reg[j])  reg_2[k] = reg[j++];
                            else reg_2[k] = cd[i++];
                        }
			for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 0 && rank != 0){
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank - 1, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);
                        int i, j, k;
                        for(i = j = k = procs_data_num - 1; k >= 0; k--){
                            if(cd[i] >= reg[j])	 reg_2[k] = cd[i--];
                            else reg_2[k] = reg[j--];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
		else{	}
            }
            delete [] reg;
            delete [] reg_2;
        }
	else{
            reg = new float [procs_data_num + 1];
            reg_2 = new float [procs_data_num];
            for (int i = 0; i < size / 2 + 1; i++) {
                //Even Phase
                if(rank % 2 == 0 && rank == n % size - 1){
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                        //Even Rank Sent Message To Odd Rank
                        MPI_Recv(reg, procs_data_num - 1, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //Even Rank Recv Message From Odd Rank
                        int i, j, k;
                        for(i = j = k = 0; k < procs_data_num; k++){
                            if(j == procs_data_num - 1) reg_2[k] = cd[i++];
                            else if(cd[i] > reg[j]) reg_2[k] = reg[j++];
                            else    reg_2[k] = cd[i++];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 0 && rank != size - 1){
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                        //Even Rank Sent Message To Odd Rank
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //Even Rank Recv Message From Odd Rank
                        int i, j, k;
                        for(i = j = k = 0; k < procs_data_num; k++){
                            if(cd[i] > reg[j])  reg_2[k] = reg[j++];
                            else    reg_2[k] = cd[i++];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 1 && rank == n % size){
                        MPI_Recv(reg, procs_data_num + 1, MPI_FLOAT, rank - 1, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);
                        int i, j, k;
                        for(i = k = procs_data_num - 1, j = procs_data_num; k >= 0; k--){
                            if(cd[i] > reg[j])  reg_2[k] = cd[i--];
                            else    reg_2[k] = reg[j--];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 1){
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank - 1, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);

                        int i, j, k;
                        for(i = j = k = procs_data_num - 1; k >= 0; k--){
                            if(cd[i] > reg[j])  reg_2[k] = cd[i--];
                            else    reg_2[k] = reg[j--];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else{
                }

                //Odd Phase
                if(rank % 2 == 1 && rank == n % size - 1){
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                        //Even Rank Sent Message To Odd Rank
                        MPI_Recv(reg, procs_data_num - 1, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //Even Rank Recv Message From Odd Rank
                        int i, j, k;
                        for(i = j = k = 0; k < procs_data_num; k++){
                            if(j == procs_data_num - 1)
                                    reg_2[k] = cd[i++];
                            else if(cd[i] > reg[j])
                                    reg_2[k] = reg[j++];
                            else    reg_2[k] = cd[i++];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 1 && rank != size - 1){
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                        //Even Rank Sent Message To Odd Rank
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        //Even Rank Recv Message From Odd Rank
                        int i, j, k;
                        for(i = j = k = 0; k < procs_data_num; k++){
                            if(cd[i] > reg[j])  reg_2[k] = reg[j++];
                            else    reg_2[k] = cd[i++];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 0 && rank == n % size){
                        MPI_Recv(reg, procs_data_num + 1, MPI_FLOAT, rank - 1, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);

                        int i, j, k;
                        for(i = k = procs_data_num - 1, j = procs_data_num; k >= 0; k--){
                            if(cd[i] > reg[j])  reg_2[k] = cd[i--];
                            else    reg_2[k] = reg[j--];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else if(rank % 2 == 0 && rank != 0){
                        MPI_Recv(reg, procs_data_num, MPI_FLOAT, rank - 1, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(cd, procs_data_num, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD);

                        int i, j, k;
                        for(i = j = k = procs_data_num - 1; k >= 0; k--){
                            if(cd[i] >= reg[j]) reg_2[k] = cd[i--];
                            else    reg_2[k] = reg[j--];
                        }
                        for(i = 0; i < procs_data_num; i++)
                            cd[i] = reg_2[i];
                }
                else{
                }
	    }
	}

	MPI_File F_Out;
        MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &F_Out);
        if(size > n){
                if(rank == 0)
                    MPI_File_write_at(F_Out, 0, cd, n, MPI_FLOAT, MPI_STATUS_IGNORE);
                else
                    MPI_File_write_at(F_Out, 0, cd, 0, MPI_FLOAT, MPI_STATUS_IGNORE);
        }
	else if(n % size == 0)
		MPI_File_write_at(F_Out, sizeof(float) * (procs_data_num) * rank, cd, procs_data_num, MPI_FLOAT, MPI_STATUS_IGNORE);
        else {
		if (rank < n % size)
                    MPI_File_write_at(F_Out, sizeof(float) * (procs_data_num) * rank, cd, procs_data_num, MPI_FLOAT, MPI_STATUS_IGNORE);
		else
                    MPI_File_write_at(F_Out, sizeof(float) * (procs_data_num) * rank + sizeof(float) * (n % size), cd, procs_data_num, MPI_FLOAT, MPI_STATUS_IGNORE);

	}
        MPI_Finalize();
	return 0;
}
