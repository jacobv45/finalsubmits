#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

extern double montecarlo(int m);

/*computes the value of pi using a hundred, ten thousand, and one million iterations of the montecarlo routine, and displays the relative error of each computation*/
int main(int argc, char** argv)
{
   
srand((unsigned)time(NULL));         //seed random numbers with timemn
printf("%f %f\n",montecarlo(100),fabs(montecarlo(100)-M_PI));

printf("%f %f\n",montecarlo(10000),fabs(montecarlo(10000)-M_PI));

printf("%f %f\n",montecarlo(1000000),fabs(montecarlo(1000000)-M_PI));
return 0;

}



double montecarlo(int m){
   int i,count=0;		    //i, for looping. count number of hits in the circle
   double x,y,z;        	    //pairs used for points in the unit square
for ( i=0; i<m; i++) {              //100 iterations of generating points, calculating
      x = (double)rand()/RAND_MAX;     //length and counting if less than one
      y = (double)rand()/RAND_MAX;
      z = x*x+y*y;
      if (z<=1) count++;
      }
   return (double)count/m*4;
}

