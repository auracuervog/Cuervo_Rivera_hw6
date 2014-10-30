#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define masa 1.69E-27
#define q 1.6E-19
#define rt 6378100
float *RK4(float h,float*posiciones,float*Velocidades);
float *accel(float *pos, float *v, float gamma);
float *vel(float *v);
int main()
{
  int i;
  char filename[100]="velocidades.dat";
  float h=0.1;
  float mint=0;
  float maxt=100;
  int N=((maxt-mint)/h);
  float *x;
  float *y;
  float *z;
  float *Vx;
  float *Vy;
  float *Vz;
  float *t;
  float *PosicionesVelocidades;
  FILE *IN;
  //Vectores Runge Kutta
  Vx=malloc(sizeof(float)*N);
  Vy=malloc(sizeof(float)*N);
  Vz=malloc(sizeof(float)*N);
  x=malloc(sizeof(float)*N);
  y=malloc(sizeof(float)*N);
  z=malloc(sizeof(float)*N);
  t=malloc(sizeof(float)*N);
  PosicionesVelocidades=malloc(sizeof(float)*7);
  
  IN=fopen(filename,"w");
   if(!IN)
    {
      printf("problems opening the file %s\n", filename);
      exit(1);
    }
  

   //Condiciones Iniciales posicion, Velocidad 
   float *VEL;
   float *POS;
   VEL=malloc(sizeof(float)*4);
   POS=malloc(sizeof(float)*3);
   VEL[0]=mint;
   VEL[1]=0;
   VEL[2]=2;
   VEL[3]=2;
   POS[0]=2*rt;
   POS[1]=0;
   POS[2]=0;
   
   for(i=1;i<N;i++)
     {
       PosicionesVelocidades=RK4(h,POS,VEL);
       
       x[i]=PosicionesVelocidades[0];
       y[i]=PosicionesVelocidades[1];
       z[i]=PosicionesVelocidades[2];
       t[i]=PosicionesVelocidades[3];
       Vx[i]=PosicionesVelocidades[4];
       Vy[i]=PosicionesVelocidades[5];
       Vz[i]=PosicionesVelocidades[6];

       POS[0]=PosicionesVelocidades[0];
       POS[1]=PosicionesVelocidades[1];
       POS[2]=PosicionesVelocidades[2];
       VEL[0]=PosicionesVelocidades[3];
       VEL[1]=PosicionesVelocidades[4];
       VEL[2]=PosicionesVelocidades[5];
       VEL[3]=PosicionesVelocidades[6];
           
 
       fprintf(IN,"%f \t %f \t %f \t \n",x[i],y[i],z[i]);
     }
   
  return 0;
}

float *RK4(float h,float*posiciones,float*Velocidades)
{
  float gamma=1;
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

  p1=malloc(sizeof(float)*3);
  p2=malloc(sizeof(float)*3);
  p3=malloc(sizeof(float)*3); 
  V1=malloc(sizeof(float)*4);
  V2=malloc(sizeof(float)*4);
  V3=malloc(sizeof(float)*4);
  k1p=malloc(sizeof(float)*3);
  k2p=malloc(sizeof(float)*3);
  k3p=malloc(sizeof(float)*3);
  k4p=malloc(sizeof(float)*3);
  k1V=malloc(sizeof(float)*3);
  k2V=malloc(sizeof(float)*3);
  k3V=malloc(sizeof(float)*3);
  k4V=malloc(sizeof(float)*3); 
  kavp=malloc(sizeof(float)*3);
  kavV=malloc(sizeof(float)*3);

  float *posactualvelactual;
  posactualvelactual=malloc(sizeof(float)*7);
  
  float *velprima1;
  velprima1=malloc(sizeof(float)*3);
  velprima1=accel(posiciones,Velocidades,gamma);
  float *posprima1;
  posprima1=malloc(sizeof(float)*3);
  posprima1=vel(Velocidades);

  k1p[0]=posprima1[0];
  k1p[1]=posprima1[1];
  k1p[2]=posprima1[2];
  k1V[0]=velprima1[0];
  k1V[1]=velprima1[1];
  k1V[2]=velprima1[2];
  // Primer Paso
  V1[0]=Velocidades[0]+(h/2);
  V1[1]=Velocidades[1]+(h/2)*k1V[0];
  V1[2]=Velocidades[2]+(h/2)*k1V[1];
  V1[3]=Velocidades[3]+(h/2)*k1V[2];

  p1[0]=posiciones[0]+(h/2)*k1p[0];
  p1[1]=posiciones[1]+(h/2)*k1p[1];
  p1[2]=posiciones[2]+(h/2)*k1p[2];


  float *velprima2;
  velprima2=malloc(sizeof(float)*3);
  velprima2=accel(p1,V1,gamma);
  float *posprima2;
  posprima2=malloc(sizeof(float)*3);
  posprima2=vel(V1);

  k2p[0]=posprima2[0];
  k2p[1]=posprima2[1];
  k2p[2]=posprima2[2];
  k2V[0]=velprima2[0];
  k2V[1]=velprima2[1];
  k2V[2]=velprima2[2];
  // Segundo Paso
  V2[0]=Velocidades[0]+(h/2);
  V2[1]=Velocidades[1]+(h/2)*k2V[0];
  V2[2]=Velocidades[2]+(h/2)*k2V[1];
  V2[3]=Velocidades[3]+(h/2)*k2V[2];

  p2[0]=posiciones[0]+(h/2)*k2p[0];
  p2[1]=posiciones[1]+(h/2)*k2p[1];
  p2[2]=posiciones[2]+(h/2)*k2p[2];


  float *velprima3;
  velprima3=malloc(sizeof(float)*3);
  velprima3=accel(p2,V2,gamma);
  float *posprima3;
  posprima3=malloc(sizeof(float)*3);
  posprima3=vel(V2);

  k3p[0]=posprima3[0];
  k3p[1]=posprima3[1];
  k3p[2]=posprima3[2];
  k3V[0]=velprima3[0];
  k3V[1]=velprima3[1];
  k3V[2]=velprima3[2];
  // Tercer Paso
  V3[0]=Velocidades[0]+(h);
  V3[1]=Velocidades[1]+(h)*k3V[0];
  V3[2]=Velocidades[2]+(h)*k3V[1];
  V3[3]=Velocidades[3]+(h)*k3V[2];

  p3[0]=posiciones[0]+(h/2)*k3p[0];
  p3[1]=posiciones[1]+(h/2)*k3p[1];
  p3[2]=posiciones[2]+(h/2)*k3p[2];


  float *velprima4;
  velprima4=malloc(sizeof(float)*3);
  velprima4=accel(p3,V3,gamma);
  float *posprima4;
  posprima4=malloc(sizeof(float)*3);
  posprima4=vel(V3);

  k4p[0]=posprima4[0];
  k4p[1]=posprima4[1];
  k4p[2]=posprima4[2];
  k4V[0]=velprima4[0];
  k4V[1]=velprima4[1];
  k4V[2]=velprima4[2];


  //Cuarto Paso
  kavp[0]=(k1p[0]+2*k2p[0]+2*k3p[0]+k4p[0])/6;
  kavp[1]=(k1p[1]+2*k2p[1]+2*k3p[1]+k4p[1])/6;
  kavp[2]=(k1p[2]+2*k2p[2]+2*k3p[2]+k4p[2])/6;

  kavV[0]=(k1V[0]+2*k2V[0]+2*k3V[0]+k4V[0])/6;
  kavV[1]=(k1V[1]+2*k2V[1]+2*k3V[1]+k4V[1])/6;
  kavV[2]=(k1V[2]+2*k2V[2]+2*k3V[2]+k4V[2])/6;


  posactualvelactual[0]=posiciones[0]+h*kavp[0];
  posactualvelactual[1]=posiciones[1]+h*kavp[1];
  posactualvelactual[2]=posiciones[2]+h*kavp[2];
  posactualvelactual[3]=Velocidades[0]+h;
  posactualvelactual[4]=Velocidades[0]+h*kavV[0];
  posactualvelactual[5]=Velocidades[1]+h*kavV[1];
  posactualvelactual[6]=Velocidades[2]+h*kavV[2];
  

  return posactualvelactual;  
  
}
float *accel(float *pos, float *v, float gamma)
{
  int i;
  float *a;
  float A,r,bo;
  int n;

  n = 2;
  a=malloc(sizeof(float)*n); 
  for(i=0;i<n;i++)
    {
      a[i]=0.0;
    }

  bo = 3E-5;
  r = sqrt(pow(pos[0],2)+pow(pos[1],2)+pow(pos[2],2)); 
  A = -(bo*pow(rt,3))/(pow(r,5)*masa*gamma);    

  a[0] = (q*A*((2*pow(pos[2],2))-pow(pos[0],2)-pow(pos[1],2))*v[2])-(3*q*A*pos[1]*pos[2]*v[3]);
  a[1] = -(q*A*((2*pow(pos[2],2))-pow(pos[0],2)-pow(pos[1],2))*v[1])+(3*q*A*pos[0]*pos[2]*v[3]);
  a[2] = (3*q*A*pos[0]*pos[2]*v[1])-(3*q*A*pos[0]*pos[2]*v[2]);
  
  return a;
}

float *vel(float *v)
{
  int i;
  float *dp;
  int n;

  n=3;
  dp=malloc(sizeof(float)*n); 
  for(i=0;i<n;i++)
    {
      dp[i]=0.0;
    }

  dp[0] = v[0];
  dp[1] = v[1];
  dp[2] = v[2];

  return dp;

}
