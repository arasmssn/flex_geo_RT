#include "ray.h"
#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

// the hidden verbose level defined in ray_util.c
extern int ray_util_verbose;

char errstr[2048];
char usage_str[]="usage: diffuser -s <seeing> (arcsec) < in_rays > out_rays";
int main(int argc,char *argv[]) {
  // diffuser takes one argument - seeing (in arcseconds)
  // with which each incoming ray is deflected irrespective of
  // wavelength or anything else. pretty unphysical..
  ray_util_verbose=0;
  float seeing=0.6;
  float kmod;
  int idum=-1;
  float sigma,arcsec;
  lin_trans lintrans;
  ray aray,x,y;
  vec cp,lat1,lat2;
  float deflection,theta;
  vec prev_u={1.0,0.0,0.0};
  int kolmogorov=0;

  arcsec=(1.0/3600.0)*atan2(1,1)/45;
  while (--argc) {
    argv++;
    switch(argv[0][0]) {
    case '-':
      switch(argv[0][1]) {
      case 'k':
	// specify kolmogorov form
	// this is 90% with 2D gaussian (sigma)
	// and 10% with 2D Gaussian (2*sigma)
	kolmogorov=1;
	break;
      case 's':
	// specify seeing.
	--argc;argv++;
	seeing=atof(argv[0]);
	break;
      default:
	sprintf(errstr,"didn't expect switch: %s\n",argv[0]);
	complain(errstr);
	break;
      }
      break;
    default:
      sprintf(errstr,"not a switch: %s\n",argv[0]);
      complain(errstr);
      break;
    }
  }

  {
    // seed the random generator.
    // measure time and invoke ran1() some number of times.
    struct timeval tv;
    struct timezone tz;
    int    i;
    gettimeofday(&tv,&tz);
    for (i=0;i<(tv.tv_sec+tv.tv_usec)%65536;i++) ran1(&idum);
    // there. initialized.
  }

  sigma=seeing*arcsec/(2*sqrt(2*log(2)));
  while (fread(&aray,sizeof(ray),1,stdin)) {
    // modify the ray using a gaussian deviate in each direction.
    // need to use the ray's original direction to calculate the two
    // other directions.
    vec *u=&aray.k;
    kmod=modulus(u);
    unitvec(u); // u will be reconstituted later after diffusing it.

    cross_prod(&prev_u,u,&cp);
    if (modulus(&cp) > 1e-6) {
      // need to regenerate the transformation matrix for proper ang. diffusing
      char spec[2048];
      vec dummy={0.0,0.0,0.0};
      float ry=-asin(u->x);
      float rx=atan2(u->y,u->z);
      sprintf(spec,"-t 0 0 0 -r %13.11f %13.11f 0",rx,ry);
      fill_transform_specs(spec,&lintrans);
      memcpy(&x.p,&dummy,sizeof(vec));
      memcpy(&y.p,&dummy,sizeof(vec));
      dummy.x=1;      dummy.y=0;      dummy.z=0;
      memcpy(&x.k,&dummy,sizeof(vec));
      dummy.x=0;      dummy.y=1;      dummy.z=0;
      memcpy(&y.k,&dummy,sizeof(vec));
      tran_ray(&x,&lintrans,1);
      tran_ray(&y,&lintrans,1);
    }

    memcpy(&lat1,&x.k,sizeof(vec));
    memcpy(&lat2,&y.k,sizeof(vec));

    cpvec(u,&prev_u);
    
    if (kolmogorov) {
      if (ran1(&idum)<0.1) {
	deflection=tan(2*sigma*sqrt(2*expdev(&idum)));
      } else {
	deflection=tan(sigma*sqrt(2*expdev(&idum)));
      }
    } else {
      deflection=tan(sigma*sqrt(2*expdev(&idum)));
    }
    theta=2*M_PI*ran1(&idum);

    scalevec(&lat1,sin(theta)*deflection);
    scalevec(&lat2,cos(theta)*deflection);
    vec_add(u,&lat1,u);
    vec_add(u,&lat2,u);
    unitvec(u);
    scalevec(u,kmod);
    // aray is ready to go..
    fwrite(&aray,sizeof(ray),1,stdout);
  }
  exit(0);
}
