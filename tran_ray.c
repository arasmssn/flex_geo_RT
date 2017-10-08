#include "ray.h"
#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <argz.h>

char  err[2048];
char  *usage_str=
"tran_rays [-d (degrees)] -t <tx> <ty> <tz> -r <rx> <ry> <rz>\n"
"                              < <input_rays> > <output_rays>\n";

int main (int argc,char *argv[]) {
  ray    aray;
  char *argz,transform_string[1024];
  size_t argz_len;

  lin_trans lintrans;

  argc--;  argv++;
  parse_transform_subargs(transform_string,&argc,&argv);
  fill_transform_specs(transform_string,&lintrans);
  while (fread(&aray,sizeof(ray),1,stdin)) {
    tran_ray(&aray,&lintrans,lintrans.invert?-1:1);
    // finished. output the ray.
    fwrite(&aray,sizeof(ray),1,stdout);
  }
  return(0);
}

