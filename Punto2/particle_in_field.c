#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv){
  double eko;
  double alpha;
  double masa_proton, q_proton;
  double c,vo,vy,vz;

  eko = atof(argv[1]); //MeV
  alpha = atof(argv[2]);
  masa_proton = pow(1.672621777*10,-27); //MeV/c**2
  q_proton = pow(1.602176487*10,-19); //C
  c = 299792458; //m/s, velocidad de la luz
  vo = (eko+(masa_proton*(pow(c,2))))/sqrt(masa_proton*eko*(pow(c,2))*(eko+(2*(pow(c,2))))); 
  vy = vo/(sen(alpha));
  vz = vo/(cos(alpha));


}
