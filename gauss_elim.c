#include<stdio.h>

/*solves the system of linear equations Ax=b using elimination*/
int main()
{
    int i,j,k,n=3;
    float A[3][3]={{-1,-5,2},{3,-3,4},{-1,3,-3}};
    float c=0.0;
    float x[3]={38,26,-27};
    
    /* constructing the augmented matrix */
    float B[3][4];
    
  for(i=0;i<n;i++){
  	for(j=0;j<n;j++){
  		B[i][j]=A[i][j];
  	}
  }
  for(i=0;i<n;i++){
  	B[i][3]=x[i];
  }
    
    for(j=0; j<n; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i=0; i<n; i++)
        {
            if(i>j)
            {
                c=B[i][j]/B[j][j];
                for(k=0; k<=n; k++)
                {
                    B[i][k]=B[i][k]-c*B[j][k];
                }
            }
        }
    }
    printf("The Upper Triangular Matrix is:\n");
    for(i=0;i<3;i++){
    	printf("%f %f %f %f\n",B[i][0],B[i][1],B[i][2],B[i][3]);
    }
    
	x[2]=B[2][3]/B[2][2];
	x[1]=(B[1][3]-B[1][2]*x[2])/B[1][1];
	x[0]=(B[0][3]-B[0][2]*x[2]-B[0][1]*x[1])/B[0][0];
	printf("\n x1=%f x2=%f x3=%f",x[0],x[1],x[2]);
    return(0);
}

