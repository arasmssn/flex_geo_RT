#include "ray.h"
#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

char  *usage_str=
"print_rays < <input_ray_set> \n";

int main (int argc,char *argv[]) {
  ray aray;
  long iray;
  if (--argc) complain("usage\n");
  iray=0;
  while (fread(&aray,sizeof(ray),1,stdin)) {
    if (!(isnan(aray.p.x) || isnan(aray.p.y) || isnan(aray.p.z))) {
      printf("# %6ld    k: %15.10lg %15.10lg %15.10lg p: %15.10lg %15.10lg %15.10lg\n",
	     iray,aray.k.x,aray.k.y,aray.k.z,aray.p.x,aray.p.y,aray.p.z);
      iray++;
    }
  }
  return(0);
}

