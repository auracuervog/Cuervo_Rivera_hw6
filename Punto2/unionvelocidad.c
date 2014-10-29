#include <stdio.h>
#include <stdlib.h>
double *accel(double *pos, double *v, double gamma);
int main()
{
  char filename[100]="velocidades.dat"
  float h=0.01;
  float mint=0;
  float maxt=100;
  int npoints=((maxt-mint)/h);
  float *x;
  float *y;
  float *z;
  float *Vx;
  float *Vy;
  float *Vz;
  float *t;
  float *Velocidades;
  float *Posiciones;
  FILE *IN;
  Vx=malloc(sizeof(float)*N);
  Vy=malloc(sizeof(float)*N);
  Vz=malloc(sizeof(float)*N);
  x=malloc(sizeof(float)*N);
  y=malloc(sizeof(float)*N);
  z=malloc(sizeof(float)*N);
  t=malloc(sizeof(float)*N);
  Velocidades=malloc(sizeof(float)*4);
  Posiciones=malloc(sizeof(float)*3);

  in=fopen(filename,"w");
   if(!in)
    {
      printf("problems opening the file %s\n", filename);
      exit(1);
    }
   for(i=1;i<npoints;i++)
     {
       Velocidades=RK4(t[i-1],x[i-1],y[i-1],z[i-1],h);
       t[i]=Velocidades[0];
       Vx[i]=Velocidades[1];
       Vy[i]=Velocidades[2];
       Vz[i]=Velocidades[3];
       fprintf(in,"%f \t %f \t %f \t %f \t \n",t[i],Vx[i],Vy[i],Vz[i]);
     }
  return 0;
}

float *RK4(float h,float*posiciones,float*velocidades)
{
  //vector de posiciones x,y,z...vector de velocidades t,vx,vy,vz
  float*p1;
  float*V1;
  float*p2;
  float*V2;
  float*p3;
  float*V3;
  float*k1p;
  float*k2p;
  float*k3p;
  float*k4p;
  float*k1V;
  float*k2V;
  float*k3V;
  float*k4V;
  float*kavp;
  float*kavV;

  float *posactual;
  float *velactual;

  float *velprima1=accel(posiciones,velocidades);
  float *posprima1=vel(velocidades);
  k1x=posprima1[0];
  k1y=posprima1[1];
  k1z=posprima1[2];
  k1vx=velprima1[0];
  k1vy=velprima1[1];
  k1vz=velprima1[2];
  // Primer Paso
  t1=velocidades[0]+(h/2);
  x1=posiciones[0]+(h/2)*k1x;
  y1=posiciones[1]+(h/2)*k1y;
  z1=posiciones[2]+(h/2)*k1z;
  Vx1=velocidades[1]+(h/2)*k1vx;
  Vy1=velocidades[2]+(h/2)*k1vy;
  Vz1=velocidades[3]+(h/2)*k1vz;
  float *prima2=accel(t1,x1,y1,z1);
  k2x=prima2[0];
  k2y=prima2[1];
  k2z=prima2[2];
  //Segundo Paso
  t2=told+(h/2);
  x2=xold+(h/2)*k2x;
  y2=yold+(h/2)*k2y;
  z2=zold+(h/2)*k2z;
  float *prima3=accel(t2,x2,y2,z2);
  k3x=prima3[0];
  k3y=prima3[1];
  k3z=prima3[2];
  //Tercer Paso
  t3=told+(h);
  x3=xold+(h)*k3x;
  y3=yold+(h)*k3y;
  z3=zold+(h)*k3z;
  float *prima4=accel(t3,x3,y3,z3);
  k4x=prima4[0];
  k4y=prima4[1];
  k4z=prima4[2];

  //Cuarto Paso
  kavx=(k1x+2*k2x+2*k3x+k4x)/6;
  kavy=(k1y+2*k2y+2*k3y+k4y)/6;
  kavx=(k1z+2*k2z+2*k3z+k4z)/6;
  
  posactual[0]=told+h;
  posactual[1]=xold+h*kavx;
  posactual[2]=yold+h*kavy;
  posactual[3]=zold+h*kavz;


  return posactual;  
  
}
double *accel(double *pos, double *v, double gamma)
{
  int i;
  double *a;
  double A,r,bo;
  int n;

  n = 2;
  a=malloc(sizeof(double)*n); 
  for(i=0;i<n;i++)
    {
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

double *vel(double *v)
{
  int i;
  double *dp;
  int n;

  n=3;
  dp=malloc(sizeof(double)*n); 
  for(i=0;i<n;i++)
    {
      dp[i]=0.0;
    }

  dp[0] = v[0];
  dp[1] = v[1];
  dp[2] = v[2];

  return dp;

}
