/*Currently produces relatively stable orbits, but the processional period is ~580 years*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define m_e 3.003467E-6 
#define m_s 1
#define G 39.44
#define I1 1.7989E-15
#define I3 1.8047E-15
#define Tday 1/365.25
#define Theta0 0.40928
#define C0 2299.933

extern double ax(double x, double y, double z, double Sx, double Sy, double Sz, double theta, double phi);
extern double ay(double x, double y, double z, double Sx, double Sy, double Sz, double theta, double phi);
extern double az(double x, double y, double z, double Sx, double Sy, double Sz, double theta, double phi);
extern double aphi(double x, double y, double z, double Sx, double Sy, double Sz, double theta,double vtheta, double phi,double vphi);
extern double atheta(double x, double y, double z, double Sx, double Sy, double Sz, double theta,double vtheta, double phi,double vphi);

int main(){
double x=1.01671388,y=0,z=0;
double Sx,Sy,Sz;
double vx=0,vy=6.177,vz=0;
double theta=0.40928, vtheta=0;
double phi=0, vphi=0;

double k1x,k2x,k3x,k4x;
double k1y,k2y,k3y,k4y;
double k1z,k2z,k3z,k4z;
double k1p,k2p,k3p,k4p;
double k1t,k2t,k3t,k4t;

double dt=0.00001;
double t=0.;
int i;

for(i=0;i<100000;i++){
Sx=-1*m_e*x/m_s;
Sy=-1*m_e*y/m_s;
Sz=-1*m_e*z/m_s;

k1x=ax(x,y,z,Sx,Sy,Sz,theta,phi);
k2x=ax(x+dt*vx/2+dt*dt*k1x/8,y,z,Sx,Sy,Sz,theta,phi);
k3x=ax(x+dt*vx/2+dt*dt*k2x/8,y,z,Sx,Sy,Sz,theta,phi);
k4x=ax(x+dt*vx+dt*dt*k3x/2,y,z,Sx,Sy,Sz,theta,phi);

k1y=ay(x,y,z,Sx,Sy,Sz,theta,phi);
k2y=ay(x,y+dt*vy/2+dt*dt*k1y/8,z,Sx,Sy,Sz,theta,phi);
k3y=ay(x,y+dt*vy/2+dt*dt*k2y/8,z,Sx,Sy,Sz,theta,phi);
k4y=ay(x,y+dt*vy+dt*dt*k3y/2,z,Sx,Sy,Sz,theta,phi);

k1z=az(x,y,z,Sx,Sy,Sz,theta,phi);
k2z=az(x,y,z+dt*vz+dt/2*dt*k1z/8,Sx,Sy,Sz,theta,phi);
k3z=az(x,y,z+dt*vz/2+dt*dt*k2z/8,Sx,Sy,Sz,theta,phi);
k4z=az(x,y,z+dt*vz+dt*dt*k3z/2,Sx,Sy,Sz,theta,phi);

k1p=aphi(x,y,z,Sx,Sy,Sx,theta,vtheta,phi,vphi);
k2p=aphi(x,y,z,Sx,Sy,Sx,theta,vtheta,phi+vphi*dt/2+dt*dt*k1p/8,vphi+dt*k1p/2);
k3p=aphi(x,y,z,Sx,Sy,Sx,theta,vtheta,phi+vphi*dt/2+dt*dt*k2p/8,vphi+dt*k2p/2);
k4p=aphi(x,y,z,Sx,Sy,Sx,theta,vtheta,phi+vphi*dt+dt*dt*k3p/2,vphi+dt*k3p);

k1t=atheta(x,y,z,Sx,Sy,Sx,theta,vtheta,phi,vphi);
k2t=atheta(x,y,z,Sx,Sy,Sx,theta+vtheta*dt/2+dt*dt*k1t/8,vtheta+dt*k1t/2,phi,vphi);
k3t=atheta(x,y,z,Sx,Sy,Sx,theta+vtheta*dt/2+dt*dt*k2t/8,vtheta+dt*k2t/2,phi,vphi);
k4t=atheta(x,y,z,Sx,Sy,Sx,theta+vtheta*dt+dt*dt*k3t/2,vtheta+dt*k3p,phi,vphi);

x+=dt*vx+dt*dt*(k1x+k2x+k3x)/6;
vx+=dt*(k1x+2*k2x+2*k3x+k4x)/6;

y+=dt*vy+dt*dt*(k1y+k2y+k3y)/6;
vy+=dt*(k1y+2*k2y+2*k3y+k4y)/6;

z+=dt*vz+dt*dt*(k1z+k2z+k3z)/6;
vz+=dt*(k1z+2*k2z+2*k3z+k4z)/6;

phi+=dt*vphi+dt*dt*(k1p+k2p+k3p)/6;
vphi+=dt*(k1p+2*k2p+2*k3p+k4p)/6;

theta+=dt*vtheta+dt*dt*(k1t+k2t+k3t)/6;
vtheta+=dt*(k1t+2*k2t+2*k3t+k4t)/6;

t+=dt;

printf("%f %f %f %f %f\n",t,x,y,phi,theta);
}

return 0;
}

double ax(double x, double y, double z, double Sx, double Sy, double Sz, double theta, double phi){
double d = sqrt((Sx-x)*(Sx-x)+(Sy-y)*(Sy-y)+(Sz-z)*(Sz-z));
double d3 = (Sx-x)*sin(theta)*sin(phi)-(Sy-y)*sin(theta)*cos(theta)+(Sz-z)*cos(theta);
double Fx = -1*(-1*G*m_s*(Sx-x)/(d*d*d)+3*G*m_s*(I1-I3)/(m_e*d*d*d*d*d)*((Sx-x)/2-5*d3*d3*(Sx-x)/(d*d)+d3*sin(phi)*sin(theta)));
return Fx;
}

double ay(double x,double y, double z, double Sx, double Sy, double Sz, double theta, double phi){
double d = sqrt((Sx-x)*(Sx-x)+(Sy-y)*(Sy-y)+(Sz-z)*(Sz-z));
double d3 = (Sx-x)*sin(theta)*sin(phi)-(Sy-y)*sin(theta)*cos(theta)+(Sz-z)*cos(theta);
double Fy = -1*(-1*G*m_s*(Sy-y)/(d*d*d)+3*G*m_s*(I1-I3)/(m_e*d*d*d*d*d)*((Sy-y)/2-5*d3*d3*(Sy-y)/(d*d)-d3*cos(phi)*sin(theta)));
return Fy;
}

double az(double x,double y,double z,double Sx,double Sy,double Sz,double theta,double phi){
double d = sqrt((Sx-x)*(Sx-x)+(Sy-y)*(Sy-y)+(Sz-z)*(Sz-z));
double d3 = (Sx-x)*sin(theta)*sin(phi)-(Sy-y)*sin(theta)*cos(theta)+(Sz-z)*cos(theta);
double Fz = -1*(-1*G*m_s*(Sz-z)/(d*d*d)+3*G*m_s*(I1-I3)/(m_e*d*d*d*d*d)*((Sz-z)/2-5*d3*d3*(Sz-z)/(d*d)+d3*cos(theta)));
return Fz;
}

double aphi(double x, double y,double z, double Sx,double Sy,double Sz, double theta, double vtheta, double phi, double vphi){
double d=sqrt((Sx-x)*(Sx-x)+(Sy-y)*(Sy-y)+(Sz-z)*(Sz-z));
double d3=(Sx-x)*sin(theta)*sin(phi)-(Sy-y)*sin(theta)*cos(theta)+(Sz-z)*cos(theta);
double Fphi=(1/sin(theta))*(3*G*m_s*d3*(I1-I3)/(I1*d*d*d*d*d)*((Sx-x)*cos(phi)+(Sy-y)*sin(phi)));
return Fphi;
}

double atheta(double x, double y,double z,double Sx, double Sy, double Sz,double theta, double vtheta, double phi, double vphi){
double d=sqrt((Sx-x)*(Sx-x)+(Sy-y)*(Sy-y)+(Sz-z)*(Sz-z));
double d3=(Sx-x)*sin(theta)*sin(phi)-(Sy-y)*sin(theta)*cos(theta)+(Sz-z)*cos(theta);
double ftheta;
ftheta = vphi*vphi*sin(theta)*cos(theta)-I3*C0*vphi*sin(theta)/I1+3*G*m_s*(I1-I3)/(I1*d*d*d*d*d)*d3*((Sx-x)*sin(phi)*cos(theta)-(Sy-y)*cos(phi)*cos(theta)-(Sz-z)*sin(theta)); 
return ftheta;
}
