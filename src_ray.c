#include "ray.h"
#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#define OA_S_X0 2500

/* based on max:min wavelengths of 29:18 */
#define OA_P_X0 1500
#define OA_P_Z0 OA_P_X0

float plnsurf (float x,float y,float z);
float offax_parabola (float x,float y,float z);
float offax_sphere (float x,float y,float z);
float local_d_spacing(float x);

long  iduml1=-1,iduml2=-1;
int   idumi;
char  *usage_str=
"ray -f <ph s-1 cm-2> -s <src x> <src y> -t <exposure time [s]>\n"
"    (-r <aperture radius [mm]> || (-x <xdim> && -y <ydim>))\n"
"    -a <cone radius [deg]> -b <lam1> <lam2> <phot_index> -- [expects stdin]\n"
"    [if -a used then -f units are ph s-1 cm-2 sr-1]\n"
"    [if -b used then -f units are ph s-1 cm-2 keV-1 (@1keV) and\n"
"         photons with wavelength lam1 < lam < lam2 are generated]\n"
"    [if -a and -b are used then the intersection of the two above cases.]\n";

int main (int argc,char *argv[]) {

  char  err[2048];
  vec   src;
  float src_flux,r_out,r_in,t_int,wavevec,xdim,ydim;
  float n_sr,area,cone_radius;
  float lam1,lam2,ph_index,np,obj_dist;
  long long int num_rays;
  int   print_surface=0;
  int   list_on_stdin=0;
  float seeing=0;
  float sigma;

  /* program to generate the ray set. 
     argument list:
     -n <number of rays> -s <source x> <source y> -a <aperture radius [mm]>
  */

  src_flux=10.0;
  int use_poisson_dist=1;
  obj_dist=0;
  r_out=0.0;
  r_in=0.0;
  xdim=ydim=0.0;
  t_int=10000.0;
  wavevec=2*M_PI/21.0;
  cone_radius=0.0;
  ph_index=-99;

  src.x=0.0;
  src.y=0.0;
  src.z=1.0;

  while (--argc) {
    argv++;
    switch (argv[0][0]) {
    case '-':
      switch(argv[0][1]) {
      case '-':
	list_on_stdin=1;
	break;
      case 'e':
	use_poisson_dist=0;
	break;
      case 'o': 
	if (argc<2) complain("argument expected after -o\n");
	--argc;	++argv;
	obj_dist=atof(argv[0]);
	break;
      case 'l': 
	if (argc<2) complain("argument expected after -l\n");
	--argc;	++argv;
	wavevec=2*M_PI/atof(argv[0]);
	break;
      case 'f': 
	if (argc<2) complain("argument expected after -f\n");
	--argc;	++argv;
	src_flux=atof(argv[0]);
	break;
      case 'S':
	// specify seeing.
	if (argc<2) complain("argument expected after -S\n");
	--argc;	++argv;
	seeing=atof(argv[0]);
	sigma=(seeing/(3600*2*sqrt(2*log(2))))*atan2(1,1)/45.0;
	fprintf(stderr,"seeing %f sigma %f\n",seeing,sigma);
	break;
      case 'a': 
	if (argc<2) complain("argument expected after -a\n");
	--argc;	++argv;
	cone_radius=atof(argv[0])*M_PI/180.0;
	n_sr=2*M_PI*(1-cos(cone_radius));
	break;
      case 'p':
	print_surface=1;
	break;  
      case 's':
	if (argc<3) complain("arguments expected after -s\n");
	--argc;	++argv;
	src.x=atof(argv[0]);
	--argc;	++argv;
	src.y=atof(argv[0]);
	if (src.x*src.x + src.y*src.y > 1.0) 
	  complain("error: -s expects vector components for the source. They"
		   " must be components of a unit vector.\n");
	break;
      case 'b':
	if (argc<4) complain("arguments expected after -b\n");
	--argc;	++argv;
	lam1=atof(argv[0]);
	--argc;	++argv;
	lam2=atof(argv[0]);
	--argc;	++argv;
	ph_index=atof(argv[0]);
	break;
      case 'x':
	if (argc<2) complain("argument expected after -x\n");
	--argc;	++argv;
	xdim=atof(argv[0]);
	break;
      case 'y':
	if (argc<2) complain("argument expected after -y\n");
	--argc;	++argv;
	ydim=atof(argv[0]);
	break;
      case 'r':
	switch (argv[0][2]) {
	case 'i':
	  if (argc<2) complain("argument expected after -ri\n");
	  --argc;	++argv;
	  r_in=atof(argv[0]);
	  break;
	case 'o':
	default:
	  if (argc<2) complain("argument expected after -r\n");
	  --argc;	++argv;
	  r_out=atof(argv[0]);
	  break;
	} 
	break;
      case 't':
	if (argc<2) complain("argument expected after -t\n");
	--argc;	++argv;
	t_int=atof(argv[0]);
	break;
      default:
	sprintf(err,"unknown switch: %s\n",argv[0]);
	complain(err);
	break;
      }
      break;
    default:
      sprintf(err,"wasn't expecting: %s\n",argv[0]);
      complain(err);
      break;
    }
  }
  
  if (((r_out==0) && (xdim==0) && (ydim==0)) || 
      (((r_out==0) && (xdim*ydim==0)) ||
       ((r_out!=0) && (xdim!=0 || ydim!=0)))) {
    complain("you need to specify the aperture. EITHER -r OR -x & -y "
	     "switches.\n");
  }
  if ((r_out!=0) && (r_in>=r_out)) {
    complain("r_in must be less than r_out!\n");
  }

  src.z=sqrt(1-src.x*src.x-src.y*src.y);
  src.x *= wavevec;
  src.y *= wavevec;
  src.z *= wavevec;

  /* cast r_out into mm. */
  {
    struct timeval tv;
    struct timezone tz;
    
    gettimeofday(&tv,&tz);
    iduml1 = tv.tv_sec + tv.tv_usec;
    iduml2 = tv.tv_sec/2 + tv.tv_usec;
    idumi=-1;
  }
  
  {
    int i,n;
    n=iduml1%4096;
    for (i=0;i<n;i++) {
      poidev(0.5,&idumi);
      ran2(&iduml1);
      ran2(&iduml2);
      ran2(&iduml2);
    }
  }

  if (r_out) {
    area=M_PI*(pow(r_out/10.0,2.0)-pow(r_in/10.0,2.0));
  } else {
    area=xdim*ydim/100.0;
  }

  if (!list_on_stdin && (cone_radius==0.0)) {
    area *= fabs(src.z/modulus(&src));
  }

  if (list_on_stdin) {
    float lam,r,phi;
    long long int i_ray;
    ray aray;
    
    //    aray.pl=0;
    //    aray.n = 1;
    while (!feof(stdin) &&
	   scanf("%g %g %lg %lg\n",&lam,&src_flux,&src.x,&src.y)) {
      fprintf(stderr,"%f %f %lf %lf\n",lam,src_flux,src.x,src.y);

      if (use_poisson_dist) {
	num_rays=floor(poidev(area*sqrt(1-(src.x*src.x+src.y*src.y))*
			      src_flux*t_int,&idumi));
      } else {
	num_rays=floor(0.5+area*sqrt(1-(src.x*src.x+src.y*src.y))*
		       src_flux*t_int);
      }

      if (num_rays) {
	wavevec=2*M_PI/lam;
	src.z=sqrt(1-src.x*src.x-src.y*src.y);
	src.x *= wavevec;
	src.y *= wavevec;
	src.z *= wavevec;
        aray.k.x=src.x;
	aray.k.y=src.y;
	aray.k.z=src.z;
	for (i_ray=0;i_ray<num_rays;i_ray++) {

	  if (i_ray%65536==65536-1) 
	    fprintf(stderr,"\r%lld/%lld (%f%% complete)..",
		    i_ray,num_rays,100*i_ray/((float)(num_rays)));
	  
	  if (r_out!=0) {
	    do {
	      r=r_out*
		sqrt((1-pow(r_in/r_out,2.0))*ran2(&iduml1)+pow(r_in/r_out,2.0));
	    } while (r==0);
	    phi=2*M_PI*ran2(&iduml2);
	    aray.p.x=r*sin(phi);
	    aray.p.y=r*cos(phi);
	    aray.p.z=0.0;
	  } else {
	    aray.p.x=xdim*(ran2(&iduml1)-0.5);
	    aray.p.y=ydim*(ran2(&iduml2)-0.5);
	    aray.p.z=0.0;
	  }
	  if (obj_dist != 0) {
	    vec p0;
	    p0.x=p0.y=0;
	    p0.z=-obj_dist;
	    cpvec(&aray.p,&aray.k);
	    vec_diff(&aray.k,&p0,&aray.k);
	    unitvec(&aray.k);
	    if (p0.z<0) {
	      scalevec(&aray.k,wavevec);
	    } else {
	      scalevec(&aray.k,-wavevec);
	    }
	  } 
	  if (seeing!=0.0) {
	    // generate lateral vector components for use later
	    
	  }
	  fprintf(stderr,"seeing is %f\n",seeing);
	  if (seeing!=0.0) {
	    // use lateral vector components and roll dice to 
	    // impart "seeing" on all rays from the source.
	    float deflection=tan(sigma*sqrt(2*expdev(&idumi)));
	    float theta=2*M_PI*ran1(&idumi);
	    vec   *avec=&aray.k;
	    float kmod=modulus(avec);
	    unitvec(avec);
	    vec   lat1={1.0,0.0,0.0};
	    vec   lat2={0.0,1.0,0.0};
	    scalevec(&lat1,deflection*sin(theta));
	    scalevec(&lat2,deflection*cos(theta));
	    vec_add(avec,&lat1,avec);
	    vec_add(avec,&lat2,avec);
	    unitvec(avec);
	    scalevec(avec,kmod);
	    //	    fprintf(stderr,"%s\n",show_vector("avec: ",avec));
	  }
	  fwrite(&aray,sizeof(ray),1,stdout);
	}
      }
    }
    exit(0);
  }

  // else continue with the regular part.

  if (ph_index != -99) {
    if (ph_index==1) {
      np=log(lam2/lam1);
    } else {
      np=1/(1-ph_index)*
	(pow(12.398/lam1,(1-ph_index))-pow(12.398/lam2,(1-ph_index)));
    }
    fprintf(stderr,"flux scaler %f\n",np);
  }

  fprintf(stderr,"source setup: %s\n",show_vector("source",&src));

  if (cone_radius>0) {
    if (ph_index!=-99) {
      if (use_poisson_dist) {
	num_rays=floor(poidev(area*src_flux*np*n_sr*t_int,&idumi));
      } else {
	num_rays=floor(0.5+area*src_flux*np*n_sr*t_int);
      }
      fprintf(stderr,"input rays: %lld between %f and %fA "
	      "(flux: %g ph s-1 cm-2 sr-2 keV-1 (@1keV); t_int %f)\n",
	      num_rays,lam1,lam2,src_flux,t_int);
    } else {
      if (use_poisson_dist) {
	num_rays=floor(poidev(area*src_flux*n_sr*t_int,&idumi));
      } else {
	num_rays=floor(0.5+area*src_flux*n_sr*t_int);
      }
      fprintf(stderr,"input rays: %lld (flux: %g ph s-1 cm-2 sr-1; t_int %f)\n",
	      num_rays,src_flux,t_int);
    }
  } else {
    if (ph_index!=-99) {
      if (use_poisson_dist) {
	num_rays=floor(poidev(area*src_flux*np*t_int,&idumi));
      } else {
	num_rays=floor(0.5+area*src_flux*np*t_int);
      }
      fprintf(stderr,"input rays: %lld (of %g expected) between %f and %fA "
	      "(flux: %g ph s-1 cm-2 keV-1 (@1keV); t_int %f)\n",
	      num_rays,
	      area*src_flux*np*t_int,
	      lam1,lam2,src_flux,t_int);
    } else {
      if (use_poisson_dist) {
	num_rays=floor(poidev(area*src_flux*t_int,&idumi));
      } else {
	num_rays=floor(0.5+area*src_flux*t_int);
      }
      fprintf(stderr,"input rays: %lld (of %g expected) (flux: %g ph s-1 cm-2; t_int %f)\n",
	      num_rays,
	      area*src_flux*t_int,
	      src_flux,t_int);
    }
  }

  if (obj_dist!=0) {
    fprintf(stderr,"obj distance : %lf\n",obj_dist);
  }

  {
    long long int i_ray;
    ray    aray;
    double r,phi;

    //    aray.pl=0;
    //    aray.n =1;
    if (cone_radius>0.0) {
      // generate unit vectors perp. to src ray
      
    }

    for (i_ray=0;i_ray<num_rays;i_ray++) {

      if (i_ray%65536==65536-1) 
	fprintf(stderr,"\r%lld/%lld (%f%% complete)..",
		i_ray,num_rays,100*i_ray/((float)(num_rays)));

      if (ph_index!=-99) {
	float cum,e1,e2,energy;
	// dial up an energy.
	cum=ran2(&iduml1);
	e1=12.398/lam2;
	e2=12.398/lam1;
	if (ph_index==1) {
	  energy=exp(cum*log(e2)+(1-cum)*log(e1));
	} else {
	  energy=exp(1/(1-ph_index)*
		     log(cum*pow(e2,(1-ph_index))+
			 (1-cum)*pow(e1,(1-ph_index))));
	}
	wavevec=2*M_PI/(12.398/energy);
      }

      /* generate the rays and dump to standard output. */
      /* draw random deviate for a circular aperture. */
      if (cone_radius>0.0) {
	// random polar angle. throw out in the limit of r=pi/2 because 
	// the projected aperture scales with cos(r)
	r=acos(1.0-n_sr*ran2(&iduml2)/(2*M_PI));
	// this selection comes in because a factor of cos(r) should
	// be used to compute survivability of a ray passing through the
	// aperture; another cos(r) because the ray is already passing
	// through the aperture. does this make sense? it gives the 
	// right answer.
	if (ran2(&iduml2) > cos(r)*cos(r)) goto next_ray;
	phi=2*M_PI*ran2(&iduml1);
	aray.k.x=sin(r)*sin(phi)*wavevec;
	aray.k.y=sin(r)*cos(phi)*wavevec;
	aray.k.z=cos(r)*wavevec;
      } else {
	unitvec(&src);
	cpvec(&src,&aray.k);
	scalevec(&aray.k,wavevec);
      }

      if (r_out!=0) {
	r=r_out*
	  sqrt((1-pow(r_in/r_out,2.0))*ran2(&iduml1)+pow(r_in/r_out,2.0));
	phi=2*M_PI*ran2(&iduml2);
	aray.p.x=r*sin(phi);
	aray.p.y=r*cos(phi);
	aray.p.z=0.0;
      } else {
	aray.p.x=xdim*(ran2(&iduml1)-0.5);
	aray.p.y=ydim*(ran2(&iduml2)-0.5);
	aray.p.z=0.0;
      }
      if (obj_dist != 0) {
	vec p0;
	p0.x=p0.y=0;
	p0.z=-obj_dist;
	cpvec(&aray.p,&aray.k);
	vec_diff(&aray.k,&p0,&aray.k);
	unitvec(&aray.k);
	if (p0.z<0) {
	  scalevec(&aray.k,wavevec);
	} else {
	  scalevec(&aray.k,-wavevec);
	}
      }
      
      if (seeing!=0.0) {
	// use lateral vector components and roll dice to 
	// impart "seeing" on all rays from the source.
	float deflection=tan(sigma*sqrt(2*expdev(&idumi)));
	float theta=2*M_PI*ran1(&idumi);
	vec   *avec=&aray.k;
	float kmod=modulus(avec);
	unitvec(avec);
	vec   lat1={1.0,0.0,0.0};
	vec   lat2={0.0,1.0,0.0};
	scalevec(&lat1,deflection*sin(theta));
	scalevec(&lat2,deflection*cos(theta));
	vec_add(avec,&lat1,avec);
	vec_add(avec,&lat2,avec);
	unitvec(avec);
	scalevec(avec,kmod);
	//	fprintf(stderr,"%s\n",show_vector("avec: ",avec));
      }
      fwrite(&aray,sizeof(ray),1,stdout);
    next_ray:
      continue;
    }
  }
  return(0);
}
