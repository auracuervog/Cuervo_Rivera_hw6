#include <stdio.h>
#include <stdlib.h>
int main()
{
  char filename[100]="velocidades.dat"
  float h=0.01;
  float mint=0;
  float maxt=100;
  int npoints=((maxt-mint)/h);
  float *Vx;
  float *Vy;
  float *Vz;
  float *t;
  float *Velocidades;
  FILE *IN;
  Vx=malloc(sizeof(float)*N);
  Vy=malloc(sizeof(float)*N);
  Vz=malloc(sizeof(float)*N);
  t=malloc(sizeof(float)*N);
  Velocidades=malloc(sizeof(float)*4);

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

float *RK4(float h,float told,float xold,float yold, float zold)
{
  float t1;
  float x1;
  float y1;
  float z1;
  float t2;
  float x2;
  float y2;
  float z2;
  float t3;
  float x3;
  float y3;
  float z3;
  float k1x;
  float k2x;
  float k3x;
  float k4x;
  float k1y;
  float k2y;
  float k3y;
  float k4y;
  float k1z;
  float k2z;
  float k3z;
  float k4z;
  float kavx;
  float kavy;
  float kavz;
  float *posactual;
  float *prima1=accel(told,xold,yold,zold);
  
  k1x=prima1[0];
  k1y=prima1[1];
  k1z=prima1[2];
  // Primer Paso
  t1=told+(h/2);
  x1=xold+(h/2)*k1x;
  y1=yold+(h/2)*k1y;
  z1=zold+(h/2)*k1z;
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
