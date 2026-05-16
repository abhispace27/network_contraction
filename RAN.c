#include <stdio.h>
#include <math.h>

/**********************
 * random # generator *
 **********************

 RANSEED(int,int) seeds the random number generator. Need to do this first.

 RAN()  returns a double between 0 and 1.

 GRAN() returns a gaussian random number with a sigma of 1.

*/
 
double c,cd,cm,u[97];
int i97,j97 ;

void RANSEED(int ij,int kl)
{
  int i,ii,j,jj,k,l,m ;
  double s,t ;

  i=((ij/177)%177)+2 ;
  j=(ij%177)+2 ;
  k=((kl/169)%178)+1 ;
  l=kl%169 ;
  for (ii=0;ii<97;ii++) {
    s=0.0 ;
    t=0.5 ;
    for (jj=0;jj<24;jj++) {
      m=(((i*j)%179)*k)%179 ;
      i=j;
      j=k;
      k=m;
      l=(53*l+1)%169;
      if (((l*m)%64)>=32) s+=t;
      t*=0.5;
    }
    u[ii]=s;
  }
  c=362436.0/16777216.0;
  cd=7654321.0/16777216.0;
  cm=16777213.0/16777216.0;
  i97=96;
  j97=32;
}
    
double rand_double()
{
  double uni;

  uni=u[i97]-u[j97];
  if (uni<0.0) uni+=1.0;
  u[i97]=uni;
  if (--i97<0) i97=96;
  if (--j97<0) j97=96;
  c-=cd;
  if (c<0.0) c+=cm;
  uni-=c;
  if (uni<0.0) uni+=1.0;
  return(uni);
}

float GRAN()
{
  float x,r;
  static float y;
  static int flag=0;
static double pi=3.141592653589793238462643;

  if(flag==0){
    x=rand_double();
    y=2*pi*rand_double();
    r=sqrt(-2.0*log(x));
    x=r*cos(y);
    y=r*sin(y);
    flag=1;
  }
  else{
    x=y;
    flag=0;
  }
  return x;
}
