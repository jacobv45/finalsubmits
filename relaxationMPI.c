#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

# define ROWS 18
# define COLS 5


/*solves the Poisson equation for a 16x5 grid using the method of relaxation*/

int main(int argc, char **argv){
int rank,size;

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&size);
MPI_Status status;

int i,j,k,m;
int n=300;      	/*number if iterations of the relaxation*/
int i_first,i_last;

double local[ROWS][COLS],local_new[ROWS][COLS];

/*procs figuring out their first and last rows or aborting if wrong # of procs*/
if(size!=2 && size!=4){
	MPI_Abort(MPI_COMM_WORLD,1);
}
else{
	m=(ROWS-2)/size;
	i_first=m*rank+1;
	i_last=m*(rank+1);
}

/*both procs fill their part of the grid(s)*/
for(i=0;i<ROWS;i++){
        for(j=0;j<COLS;j++){
                local_new[i][j]=0;
                local[i][j]=0;
        }
}
/*and dont forget our source point!*/
if(rank==0){
local[0][2]=1;
local_new[0][2]=1;
}


/*the relaxation part ONLY ON THE INTERIOR POINTS*/

for(k=0;k<n;k++){
        for(i=i_first;i<=i_last;i++){
                for(j=1;j<COLS-1;j++){
                        local_new[i][j]=(local[i][j-1]+local[i][j+1]+local[i-1][j]+local[i+1][j])/4;
                }
        }
memmove(*local,*local_new,COLS*ROWS*sizeof(double));
/*send down if not on bottom*/	
	if(rank<size-1){
	MPI_Send(local[i_last],COLS,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
	MPI_Recv(local[i_last+1],COLS,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,&status);
	}	
/*send up if not on top */
	if(rank>0){
	MPI_Send(local[i_first],COLS,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);
	MPI_Recv(local[i_first-1],COLS,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status);
	}
MPI_Barrier(MPI_COMM_WORLD);
}

/*print out the results*/
for(k=0;k<size;k++){
if(rank==k){
	for(i=i_first;i<=i_last;i++){
		for(j=0;j<COLS;j++){
			printf("%f\t",local[i][j]);
		}	
		puts("");
	}
}
}


MPI_Finalize();
return 0;
}

