#include <stdio.h>
#include <math.h>
#include <stdlib.h>

extern double up(double x,int n);

int main(void){

    double x;
    int n;

    printf("Enter x and n:");
    scanf("%lf%d", &x, &n);

    printf("x = %f\t j_%d = %f",x,n,up(x,n));
return(0);
}

double up(double x, int n){
	double one, two, thr;
	int k;

	one=(sin(x))/x;
	two=(sin(x)-x*cos(x))/(x*x);

for(k=1;k<n;k++){
	thr=((2.*k+1.)/x)*two-one;
	one=two;
	two=thr;
}
return(thr);
}
