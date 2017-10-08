#include <stdlib.h>
#include <stdio.h>
#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>

#define ZERNIKE_HOME
#include "zernike.h"
#undef ZERNIKE_HOME

double zernike (float rho,float theta,int n,int m,
		float *dz_drho,float *dz_dtheta) {
  if (n<0) {
    fprintf(stderr,"error - zernike index n not legal: n=%d\n",n);
    exit(1);
  }
  if (abs(m)>n || ((n-abs(m))%2!=0)) {
    fprintf(stderr,"error - zernike index m for given n(%d) not legal: m=%d\n",n,m);
    exit(1);
  }
  {
    int s,s_max=(n-abs(m))/2;
    double Rnm=0,Nnm=0,dRnm_drho=0;
    float amp;
    if (z_edge_amplitudes==0) {
      Nnm=sqrt(2*(n+1)/(m==0?2:1));
    } else {
      Nnm=1;
    }
    for (s=0;s<=s_max;s++) {
      amp=(pow(-1,s)*factrl(n-s))/
	(factrl(s)*factrl(abs(m)+s_max-s)*factrl(s_max-s));
      Rnm+=amp*pow(rho,n-2*s);
      if (n-2*s>0) {
	dRnm_drho+=amp*(n-2*s)*pow(rho,n-2*s-1);
      }
    }

    if (m>=0) {
      *dz_drho=Nnm*dRnm_drho*cos(m*theta);
      *dz_dtheta=(rho==0)?0:-Nnm*Rnm*m*sin(m*theta);
      return(Nnm*Rnm*cos(m*theta));
    } else {
      *dz_drho=-Nnm*dRnm_drho*sin(m*theta);
      *dz_dtheta=(rho==0)?0:-Nnm*Rnm*m*cos(m*theta);
      return(-Nnm*Rnm*sin(m*theta));
    }
  }
}

double zernike_distortions (float rho,float theta,float *coeffs,int nz,
			    float *dz_drho_tot,float *dz_dtheta_tot) {
  int j,m,n;
  double value=0;
  float dz_drho,dz_dtheta;

  *dz_drho_tot  =0.0;
  *dz_dtheta_tot=0.0;

  for (j=0;j<nz;j++) {
    if (coeffs[j]!=0) {
      n=ceil((-3+sqrt(9+8*j))/2);
      m=2*j-n*(n+2);
      value+=coeffs[j]*zernike(rho,theta,n,m,&dz_drho,&dz_dtheta);
      *dz_drho_tot   += coeffs[j]*dz_drho;
      *dz_dtheta_tot += coeffs[j]*dz_dtheta;
    }
  }
  //  fprintf(stderr,"zern: rho %f theta %f zdist %g\n",rho,theta,value);
  return(value);
}

