#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXPNT 100 /* maximum number of points */
#define pi 3.14159
#define q 1.602176487e-19 //C,carga del proton
#define masa 1.672621777e-27 //MeV/c**2, masa del proton
#define c 299792458 //m/s, velocidad de la luz
#define rt 6378100 //m

double *accel(double *pos, double *v, double gamma);
double *vel(double *v);

int main(int argc, char **argv){
  double eko, alpha, gamma, vo;
  int i,n;  
  double *pos;
  double *v;
 
  n = 2;

  pos=malloc(sizeof(double)*n); 
  v=malloc(sizeof(double)*n);  

  /*Inicializando los arreglos*/
  for(i=0;i<n;i++){
    pos[i]=0.0;
    v[i]=0.0;
  }

  /*Condiciones iniciales:*/
  eko = atof(argv[1]); //MeV
  alpha = atof(argv[2]);

  vo = (eko+(masa*(pow(c,2))))/sqrt(masa*eko*(pow(c,2))*(eko+(2*(pow(c,2))))); 
  gamma = 1/(sqrt(1-((pow(vo,2))/(pow(c,2))))); 

  //Posicion inicial
  pos[0] = 2*rt; 
  pos[1] = 0.0;
  pos[2] = 0.0;
 
  //Velocidad inicial
  v[0] = 0.0;
  v[1] = vo*sin(alpha*(pi/180));
  v[2] = vo*cos(alpha*(pi/180));

}

double *accel(double *pos, double *v, double gamma){
  int i;
  double *a;
  double A,r,bo;
  int n;

  n = 2;
  a=malloc(sizeof(double)*n); 
  for(i=0;i<n;i++){
    a[i]=0.0;
  }

  bo = pow(3*10,-5);
  r = sqrt(pow(pos[0],2)+pow(pos[1],2)+pow(pos[2],2)); 
  A = -(bo*pow(rt,3))/(pow(r,5)*masa*gamma);    

  a[0] = (q*A*((2*pow(pos[2],2))-pow(pos[0],2)-pow(pos[1],2))*v[1])-(3*q*A*pos[1]*pos[2]*v[2]);
  a[1] = -(q*A*((2*pow(pos[2],2))-pow(pos[0],2)-pow(pos[1],2))*v[0])+(3*q*A*pos[0]*pos[2]*v[2]);
  a[2] = (3*q*A*pos[0]*pos[2]*v[0])-(3*q*A*pos[0]*pos[2]*v[1]);
  
  return a;
}
double *vel(double *v){
  int i;
  double *dp;
  int n;

  n=2;
  dp=malloc(sizeof(double)*n); 
  for(i=0;i<n;i++){
    dp[i]=0.0;
  }

  dp[0] = v[0];
  dp[1] = v[1];
  dp[2] = v[2];

  return dp;

}

