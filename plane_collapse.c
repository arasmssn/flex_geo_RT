#include "ray.h"
#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

double plnsurf (vec *v);

float xco,yco,zco,cco;
char  *usage_str=
"plane_collapse -c <xco> <yco> <zco> <cco>\n"
"               for the plane with (xco*x + yco*y + zco*z == cco)\n";

int main (int argc,char *argv[]) {

  char  err[2048];

  xco=yco=0.0;
  zco=1;
  cco=0.0;
  
  while (--argc) {
    argv++;
    switch (argv[0][0]) {
    case '-':
      switch(argv[0][1]) {
      case 'c': 
	if (argc<5) complain("arguments expected after -c\n");
	--argc;	++argv;	xco=atof(argv[0]);
	--argc;	++argv;	yco=atof(argv[0]);
	--argc;	++argv;	zco=atof(argv[0]);
	--argc;	++argv;	cco=atof(argv[0]);
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
  
  {
    ray aray;

    while (fread(&aray,sizeof(ray),1,stdin)) {
      if (ray_surface_collapse(&aray,plnsurf)) 
	fwrite(&aray,sizeof(ray),1,stdout);
    }
  }
  return(0);
}

double plnsurf (vec *v) {
  /* for now use 1*X + 0*Y + 1*Z */
  return(xco*v->x + yco*v->y + zco*v->z - cco);
}

