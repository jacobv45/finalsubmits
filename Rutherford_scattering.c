#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


int main(void){
/*beam energy in eV's, atomic radius in cm, density in kg/cm^3, atmoic number*/
double Ec=100;

double Z=79;
//double R=1.44*pow(10,-8);
double rho=0.01930;
double t=0.001;			 //target thickness in cm (0.01 would be 100 micrometers, an average metallic thin film width)

double alpha,sigma,lambda;	 //parameters needed to be calculated
//double A=M_PI*R*R;	 	 //Nucelous area
double Ei,Es,Ef;		 //Initial Energy, energy at the scattering point, final (exit) energy
double Si,Sf; 		         //Path length before scattering, path length from scattering to exit
double J;			 //mean ionization potential
double Psi,Phi;    		 //Scattering angles (psi is around Z, beam axis; phi is angle away from Z)
double z;       	         //used for smearing of incident energy
double k=0.001;          	 //used to scale z^
double dEdS;			 //rate of energy losss


srand((unsigned)time(NULL)); 	 //seeding rand with time

alpha=(0.0034)*pow(Z,0.67)/Ec;
sigma=5.21*pow(10,-24)*(Z*Z/Ec*Ec)*(4*M_PI/(alpha+alpha*alpha))*((Ec+511)/(Ec+1024))*((Ec+511)/(Ec+1024));
lambda=196/(6.02*pow(10,23)*rho*sigma);

z=sqrt(-2*log((double)rand()/RAND_MAX))*cos(2*M_PI*((double)rand()/RAND_MAX));

Ei=Ec*(1+k*z);

Si=-1*lambda*log((double)rand()/RAND_MAX);   

J=(9.76*Z+58.5/(pow(Z,0.19)))*(0.001);

dEdS=-78500*(rho*Z)*log(1+1.166*Ec/J)/(Ec*196)*1000;

Es=Ei+Si*dEdS;

double RND=(double)rand()/RAND_MAX;
double b=(2*alpha*RND)/(1+alpha-RND);

Phi=acos(1-b);

Psi=2*M_PI*((double)rand()/RAND_MAX);
/*finding Sf based on forward or backward scattering*/
if(cos(Phi)>0){
Sf=(t-Si)/cos(Phi);
}
if(cos(Phi)<0){
Sf=-Si/cos(Phi);
}

Ef=Es+Sf*dEdS;

//print results
printf("%f Incident energy (eV)\n",Ei);
printf("%f Scattering energy (eV)\n",Es);
printf("%f Exit Energy (eV)\n",Ef);
printf("%f Interaction point (nm)\n",Si*10000000);
printf("%f Distance after interaction (cm)\n",Sf);
printf("%f Phi\n",Phi);
printf("%f Psi\n",Psi);
//printf("%f %f %f\n",alpha,sigma*pow(10,16),lambda*pow(10,5));
//printf("%f %f\n",b,t-Si);
return 0;
}
