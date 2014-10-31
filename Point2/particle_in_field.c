#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define m 1.69E-27 //kg
#define q 1.60217E-19 //C
#define rt 6378100.0 //m
#define pi 3.14159
#define c 299792458 //m/s
#define bo 3E-5
float Vxprim(float Vy,float Vz);
float Vyprim(float Vx,float Vy,float Vz);
float Vzprim(float Vx,float Vy,float Vz);
float accelx(float xold, float yold, float zold, float x1old, float y1old, float z1old, float gamma, float r);
float accely(float xold, float yold, float zold, float x1old, float y1old, float z1old, float gamma, float r);
float accelz(float xold, float yold, float zold, float x1old, float y1old, float z1old, float gamma, float r);
float *RK2(float told,float xold,float yold,float zold,float Vxold,float Vyold,float Vzold,float h);

int main(int argc, char **argv)
{
  int i;
  float ek = atof(argv[1]);
  float alpha_o = atof(argv[2]);
  float *Vx;
  float *Vy;
  float *Vz;
  float *x;
  float *y;
  float *z;
  float *t;
  float *espacioyvel;
  float tmin=0;
  float tmax=20;
  float k=ek*q;
  float gamma=(k/(m*pow(c,2)))+1;
  float h=0.0001;
  FILE *in;
  float alpha=alpha_o*pi/180;
  long int N=(tmax-tmin)/h;
  espacioyvel=malloc(sizeof(float)*7);
  Vx=malloc(sizeof(float)*N);
  Vy=malloc(sizeof(float)*N);
  Vz=malloc(sizeof(float)*N);
  x=malloc(sizeof(float)*N);
  y=malloc(sizeof(float)*N);
  z=malloc(sizeof(float)*N);
  t=malloc(sizeof(float)*N);
 
  char filename[100];
  snprintf(filename, sizeof(char)*100,"trayectoria_%.1f_%.1f.dat",(float) ek, (float) alpha_o);

  float V0=c*sqrt(1-1/(pow(gamma,2)));
  t[0]=tmin;
  Vx[0]=0;
  Vy[0]=V0*sin(alpha);
  Vz[0]=V0*cos(alpha);
  x[0]=2*rt;
  y[0]=0;
  z[0]=0;

  in = fopen(filename,"w");
  if(!in)
    {
      printf("problems opening the file %s\n", filename);
      exit(1);
    }

  for(i=1;i<N;i++)
    {
      espacioyvel=RK2(t[i-1],x[i-1],y[i-1],z[i-1],Vx[i-1],Vy[i-1],Vz[i-1],h);
      t[i]=espacioyvel[0];
      Vx[i]=espacioyvel[1];
      Vy[i]=espacioyvel[2];
      Vz[i]=espacioyvel[3];
      x[i]=espacioyvel[4];
      y[i]=espacioyvel[5];
      z[i]=espacioyvel[6];
      fprintf(in,"%f \t %f \t %f \t \n",x[i],y[i],z[i]);
      //printf("%f \t %f \t %f \t %f \t \n",t[i],x[i],y[i],z[i]);
    }
  return 0;
   
}

float Vxprim(float Vx,float Vy)
{
  float Vxd;
  float sigma=10.0;
  Vxd=Vy;
  return Vxd;
}

float Vyprim(float Vx,float Vy,float Vz)
{
  float Vyd;
  float rho=28.0;
  Vyd=-Vx;
  return Vyd;
}

float Vzprim(float Vx,float Vy, float Vz)
{
  float Vzd;
  float beta=8/3;
  Vzd=0;
  return Vzd;
}

float *RK2(float told,float xold,float yold,float zold,float Vxold,float Vyold,float Vzold,float h)
{
  float Vxp1;
  float Vyp1;
  float Vzp1;
  float tmi;
  float Vxmi;
  float Vymi;
  float Vzmi;
  float Vxp2;
  float Vyp2;
  float Vzp2;
  float *posvel;
  float tnew;
  float Vxnew;
  float Vynew;
  float Vznew;
  
  float r;
  float gamma=1;

  float xp1;
  float yp1;
  float zp1;
  float xmi;
  float ymi;
  float zmi;
  float xp2;
  float yp2;
  float zp2;
  float xnew;
  float ynew;
  float znew;
  posvel=malloc(sizeof(float)*7);
  r=sqrt(pow(xold,2)+pow(yold,2)+pow(zold,2));
  Vxp1=accelx(xold,yold,zold,Vxold,Vyold,Vzold,gamma,r);
  Vyp1=accely(xold,yold,zold,Vxold,Vyold,Vzold,gamma,r);
  Vzp1=accelz(xold,yold,zold,Vxold,Vyold,Vzold,gamma,r);

  tmi=told+h/2;
  Vxmi=Vxold+(h/2)*Vxp1;
  Vymi=Vyold+(h/2)*Vyp1;
  Vzmi=Vzold+(h/2)*Vzp1;
  xmi=xold+(h/2)*Vxold;
  ymi=yold+(h/2)*Vyold;
  zmi=zold+(h/2)*Vzold;
  
  r=sqrt(pow(xmi,2)+pow(ymi,2)+pow(zmi,2));

  Vxp2=accelx(xmi,ymi,zmi,Vxmi,Vymi,Vzmi,gamma,r);
  Vyp2=accely(xmi,ymi,zmi,Vxmi,Vymi,Vzmi,gamma,r);
  Vzp2=accelz(xmi,ymi,zmi,Vxmi,Vymi,Vzmi,gamma,r);
  
  tnew=told+h;
  Vxnew=Vxold+h*Vxp2;
  Vynew=Vyold+h*Vyp2;
  Vznew=Vzold+h*Vzp2;
  xnew=xold+(h)*Vxmi;
  ynew=yold+(h)*Vymi;
  znew=zold+(h)*Vzmi;

  posvel[0]=tnew;
  posvel[1]=Vxnew;
  posvel[2]=Vynew;
  posvel[3]=Vznew;
  posvel[4]=xnew;
  posvel[5]=ynew;
  posvel[6]=znew;
  return posvel;

  
}

float accelx (float xold, float yold, float zold, float x1old, float y1old, float z1old, float gamma, float r ){
  float A;
  float x2;
  A = (bo*pow(rt,3))/(pow(r,5));    
  x2 = -(A*q/(gamma*m))*(y1old*(2*pow(zold,2)-pow(xold,2)-pow(yold,2))-3*z1old*(yold*zold));
  //x2=y1old;
  return x2;
}  
float accely (float xold, float yold, float zold, float x1old, float y1old, float z1old, float gamma, float r ){
    float A;
    float y2;
    A = (bo*pow(rt,3))/(pow(r,5));

    y2=(A*q/(gamma*m))*(x1old*(2*pow(zold,2)-pow(xold,2)-pow(yold,2))-3*z1old*(xold*zold));
    //y2=-x1old;
    return y2;
 }
float accelz (float xold, float yold, float zold, float x1old, float y1old, float z1old, float gamma, float r ){
  float A;
  float z2;
  A = (bo*pow(rt,3))/(pow(r,5));
  
  z2 = -(A*q/(gamma*m))*(x1old*(3*yold*zold)-y1old*(3*xold*zold));
  //z2=0;
  return z2;
}



