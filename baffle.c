#include "ray.h"
#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

double asp_surf (vec *v);
float n_air (float lamda);
float n_silica (float lamda);

char *usage_str="baffle -ro <r_outer> -ri <r_inner> [-P (print only)]}\n";

typedef struct asphere {
  double C;
  double k;
  float  p_atm;
  double A[11];
} asp;

typedef struct sp {
  int use;
  int arms;
  float width;
  float length;
} sp;

asp *as;
sp  spider={0,0,0.0,0.0};

int main(int argc,char *argv[]) {
  int ind,print_only=0,lens=0;
  double r_outer=4300,r_inner=1700;
  int max_transforms=10;
  int n_transform=0;
  char errstr[512],transform_strings[max_transforms][1024];
  lin_trans tf[max_transforms];
  asp asph[2];
  int obscuring_ring=0;

  for (ind=0;ind<11;ind++) {
    asph[0].A[ind]=asph[1].A[ind]=0.0;
  }

  asph[0].p_atm=asph[1].p_atm=1.0;

  asph[0].C=1/1000.0;
  asph[0].k=asph[1].k=0;

  as=&asph[0];
  as->C=0;
  as->k=0;

  while (--argc) {
    argv++;
    switch (argv[0][0]) {
    case '-':
      switch(argv[0][1]) {
      case 'T':
	if (n_transform==max_transforms) 
	  complain("too many transforms. increase max_transforms in asphere.c");
	parse_transform_subargs(transform_strings[n_transform++],&argc,&argv);
	break;
      case 'P':
	print_only=1;
	break;
      case 'r':      case 'z':      case 'Z':      case 's':      case 'A':
        if (argc<2) {
	  sprintf(errstr,"arguments expected after -%c\n",argv[0][1]);
	  complain(errstr);
	}
	switch(argv[0][1]) {
	case 's': // spider specification
	  spider.use=1;
	  switch(argv[0][2]) {
	  case 'n':
	    --argc; ++argv;
	    spider.arms=atoi(argv[0]);
	    break;
	  case 'w':
	    --argc; ++argv;
	    spider.width=atof(argv[0]);
	    break;
	  case 'l':
	    --argc; ++argv;
	    spider.length=atof(argv[0]);
	    break;
	  }
	  break;
	case 'r':
	  switch(argv[0][2]) {
	  case 'g':
	    --argc; ++argv;
	    obscuring_ring=atoi(argv[0]);
	    break;
	  case 'o':
	    --argc; ++argv; 	    
	    r_outer=atof(argv[0]);
	    break;
	  case 'i':
	    --argc; ++argv;
	    r_inner=atof(argv[0]);
	    break;
	  default:
	    complain("");
	    break;
	  }
	  break;
	case 'z':	case 'Z':
	  --argc; ++argv; 
	  (as->A)[0]=atof(argv[0]);
	  break;
	case 'A':
	  ind=atoi(argv[0]+2);
	  if (ind<1 || ind>10) 
	    complain("");
	  --argc;++argv;
	  (as->A)[ind]=atof(argv[0]);
	  (as->A)[ind] *= 1e3/pow(1e3,ind);
	  break;
	default:
	  sprintf(errstr,"this shouldn't happen\n");
	  complain(errstr);
	  break;
	}
	break;
      default:
	sprintf(errstr,"unexpected switch: %s\n",argv[0]);
	complain(errstr);
	break;
      }
      break;
    default:
      sprintf(errstr,"wasn't expecting: %s\n",argv[0]);
      complain(errstr);
      break;
    }
  }

  {
    int i=0;
    while (i<n_transform) {
      fprintf(stderr,"\nwill transform (%d) according to %s\n\n",i,transform_strings[i]);
      // here convert the transform string according to rules in tran_ray
      fill_transform_specs(transform_strings[i],&tf[i]);
      i++;
    }
  }

  as=&asph[0];

  {
    char *s;

    if (lens) {
      s="(1)";
    } else {
      s="   ";
    }
  
    fprintf(stderr,
	    "baffle setup:  \n"
	    "obsc_ring: %d\n"
	    "%s     = %+g mm\n"
	    "%s     = %+g mm\n"
	    "%s     = %+g mm\n",
	    obscuring_ring,
	    "r_outer",r_outer,
	    "r_inner",r_inner,
	    "z_plane",(as->A)[0]);

    for (ind=0;ind<11;ind++) {
      if ((as->A)[ind]!=0) {
        fprintf(stderr,"%s           A[%02d] = %+g mm^{%d}\n",
		s,ind,(as->A)[ind],1-ind);
      }
    }
    if (spider.use) {
      fprintf(stderr,
	      "spider setup:  \n"
	      "%s     = %d arms\n"
	      "%s     = %f mm\n"
	      "%s     = %f mm\n",
	      "n_arms ",spider.arms,
	      "width  ",spider.width,
	      "length ",spider.length);
    }
  }

  if (print_only) {
    vec v;
    v.z=v.y=0;
    
    as=&asph[0];
    for (v.x=-r_outer;v.x<=r_outer;v.x+=r_outer/1000) {
      if (fabs(v.x)>r_inner)
	printf("%g %g\n",v.x,asp_surf(&v));
    }
    exit(0);
  }
  {
    ray aray;
    vec inray_p;
    float r;

    as=&asph[0];
    while (fread(&aray,sizeof(ray),1,stdin)) {

      if (n_transform>0) {
	int i=0;
	while (i<n_transform) {
	  tran_ray(&aray,&tf[i],(tf[i].invert==0)?1:-1);
	  i++;
	}
      }
      
      ray_surface_collapse(&aray,asp_surf);
      memcpy(&inray_p,&aray.p,sizeof(vec));

      r=sqrt(pow(aray.p.y,2.0)+pow(aray.p.x,2.0));

      switch ((r-r_outer)*(r-r_inner)<0) {
      case 0:
	// r is outside r_inner & r_outer
	if (obscuring_ring==0) goto skip_ray;
	break;
      case 1:
	// r is between r_inner & r_outer
	if (obscuring_ring==1) goto skip_ray;
	break;
      }

      if (spider.use) {
	int arm_index;
	float x,y,theta;
	for (arm_index=0;arm_index<spider.arms;arm_index++) {
	  theta=2*M_PI*(arm_index+0.5)/(1.0*spider.arms);
	  y=aray.p.y*cos(theta)-aray.p.x*sin(theta);
	  x=aray.p.x*cos(theta)+aray.p.y*sin(theta);
	  if ((x>0) && (fabs(y)<0.5*spider.width) && 
	      ((spider.length==0) || 
	       (x<0.5*spider.length))) 
	    goto skip_ray;
	}
      }

      if (n_transform>0) {
	int i=n_transform;
	while (i>0) {
	  i--;
	  tran_ray(&aray,&tf[i],(tf[i].invert==0)?-1:1);
	}
      }

      fwrite(&aray,sizeof(ray),1,stdout);
    skip_ray:
      continue;
    }  
  }
  exit(0);
}
  
double asp_surf (vec *v) {
  double r=sqrt(pow(v->x,2.0)+pow(v->y,2.0));
  double z=as->C*pow(r,2.0)/(1+sqrt(1-(1+as->k)*pow(as->C*r,2.0)));
  int    ind;
  for (ind=0;ind<11;ind++) 
    z+=(as->A)[ind]*pow(r,ind);
  return(z-v->z);
}

