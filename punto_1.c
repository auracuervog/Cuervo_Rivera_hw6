#include <stdio.h>
#include <stdlib.h>
float xprim(float x,float y);
float yprim(float x,float y);
float *RK2(float told,float xold,float yold,float h);

int main(int argc, char **argv)
{
  char filename[100]="presa.dat";
  int i;
  float *x;
  float *y;
  float *t;
  float *espacio;
  float tmin=0;
  float tmax=1;
  float h=0.001;
  FILE *in;
  int N=(tmax-tmin)/h;
  espacio=malloc(sizeof(float)*3);
  x=malloc(sizeof(float)*N);
  y=malloc(sizeof(float)*N);
  t=malloc(sizeof(float)*N);
  float x0=atof(argv[1]);
  float y0=atof(argv[2]);

  t[0]=tmin;
  x[0]=x0;
  y[0]=y0;
  in = fopen(filename,"w");
  if(!in)
    {
      printf("problems opening the file %s\n", filename);
      exit(1);
    }

  for(i=1;i<N;i++)
    {
      espacio=RK2(t[i-1],x[i-1],y[i-1],h);
      t[i]=espacio[0];
      x[i]=espacio[1];
      y[i]=espacio[2];
      fprintf(in,"%f \t %f \t %f \t \n",t[i],x[i],y[i]);
    }
   
}

float xprim(float x,float y)
{
  float xd;
  xd=20*x-x*y;
  return xd;
}

float yprim(float x,float y)
{
  float yd;
  yd=-30*y+x*y;
  return yd;
}

float *RK2(float told,float xold,float yold,float h)
{
  float xp1;
  float yp1;
  float tmi;
  float xmi;
  float ymi;
  float xp2;
  float yp2;
  float *pos;
  float tnew;
  float xnew;
  float ynew;
  pos=malloc(sizeof(float)*2);
  xp1=xprim(xold,yold);
  yp1=yprim(xold,yold);

  tmi=told+h/2;
  xmi=xold+(h/2)*xp1;
  ymi=yold+(h/2)*yp1;
  
  xp2=xprim(xmi,ymi);
  yp2=yprim(xmi,ymi);
  
  tnew=told+h;
  xnew=xold+h*xp2;
  ynew=yold+h*yp2;

  pos[0]=tnew;
  pos[1]=xnew;
  pos[2]=ynew;
  
  return pos;

  
}
