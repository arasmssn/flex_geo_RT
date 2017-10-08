#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "numrec/include/nr.h"
// #include "numrec/include/nrutil.h"
#include "cfitsio/fitsio.h"
#include "cfitsio/fitsio2.h"
#include "ray.h"

char  err[2048];
char  *usage_str=
"image -F <flen> -l <length_scale> -d <image_dim> -o <outfile> -e [erase image first] -p <mu_per_photon> [each ray triggers a poisson diceroll and increments the pixel by this amount] -k [scale contributions by k modulus] < <input_rays> \n";

// program to produce a fits image on the fly and write it to the stdout.

void printerror( int status);

int main (int argc,char *argv[]) {

  char *outfile;
  float side,flen=1;
  int poisson_weight=1;
  int dim;
  int erase=0,invert=0;
  int kscale=0;
  int sqroot=0;
  int add_sky_background=0;
  int add_readout_noise=0;
  float poisson_mu=0;
  float readnoise_per_pixel=0;

  if (argc==1) complain("");

  while (--argc) {
    ++argv;
    switch(argv[0][0]) {
    case '-':
      switch(argv[0][1]) {
      case 's': // sky background, electrons per pixel
	--argc;argv++;
	add_sky_background=1;
	poisson_mu=atof(argv[0]);
	fprintf(stderr,"poisson_mu=%f\n",poisson_mu);
	break;
      case 'n':
	--argc;argv++;
	add_readout_noise=1;
	readnoise_per_pixel=atof(argv[0]);
	fprintf(stderr,"RON=%f\n",readnoise_per_pixel);
	break;
      case 'p': // each each photon in the bundle have this poisson weight.
	--argc;argv++;
	poisson_weight=atof(argv[0]);
	break;
      case 'o': // filename (output)
	--argc;argv++;
	outfile=argv[0];
	break;
      case 'F': // focal length, sets platescale (arcsec per pixel)
	--argc;argv++;
	flen=atof(argv[0]);
	break;
      case 'd': // image dimension (npix on a side)
	--argc;argv++;
	dim=atoi(argv[0]);
	break;
      case 'k': // scale each photon according to k vector modulus
	kscale=1;
	break;
      case 'q': // apply square root to output image
	sqroot=1;
	break;
      case 'i': // invert values (x->-x;y->-y;z->-z)
	invert=1;
	break;
      case 'e': // erase image if it exists.
	erase=1;
	break;
      case 'l': // length dimension (mm on a side)
	--argc;argv++;
	side=atof(argv[0]);
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
    fitsfile *ff;
    long naxis=2,naxes[]={dim,dim};
    long index;
    int bitpix=FLOAT_IMG,newfile,status=0;
    float *map;
    float crpix[2],crval[2],cdelt[2];
    char  *ctype[2],vals[2][512];
    int   nfound;

    if (erase) remove(outfile);
    
    if (fits_open_file(&ff,outfile,1,&status)) {
      newfile=1;
      status=0;
    } else {
      newfile=0;
    }

    if (newfile) {

      if (fits_create_file(&ff,outfile,&status))
	printerror(status);


      if (fits_create_img(ff,bitpix,naxis,naxes,&status))
	printerror(status);

	crpix[0]=dim/2;
	crval[0]=0;
	cdelt[0]=-1.0*side/(dim*flen)*180/M_PI;
	cdelt[0]=-1.0*side/(dim*flen);
	crpix[1]=dim/2;
	crval[1]=0;
	cdelt[1]=1.0*side/(dim*flen)*180/M_PI;
	cdelt[1]=1.0*side/(dim*flen);

	ctype[0]=vals[0];
	ctype[1]=vals[1];

	sprintf(ctype[0],"RA---TAN");
	sprintf(ctype[1],"DEC--TAN");

	fits_write_keys_flt(ff,"CRPIX",1,2,crpix,2,NULL,&status);
	fits_write_keys_flt(ff,"CRVAL",1,2,crval,2,NULL,&status);
	fits_write_keys_flt(ff,"CDELT",1,2,cdelt,2,NULL,&status);
	fits_write_keys_str(ff,"CTYPE",1,2,ctype,NULL,&status);
	if (status) printerror(status);
      fprintf(stderr,"this is a new file..\n");
    } else {
      fits_read_keys_lng(ff,"NAXIS",1,2,naxes,&nfound,&status);
      if (nfound!=2) complain("didn't find 2! (a)");
      fits_read_keys_flt(ff,"CRPIX",1,2,crpix,&nfound,&status);
      if (nfound!=2) complain("didn't find 2! (b)");
      fits_read_keys_flt(ff,"CRVAL",1,2,crval,&nfound,&status);
      if (nfound!=2) complain("didn't find 2! (c)");
      fits_read_keys_flt(ff,"CDELT",1,2,cdelt,&nfound,&status);
      if (nfound!=2) complain("didn't find 2! (d)");
      dim=naxes[1];
      side=flen*naxes[1]*cdelt[1];
      if (status) printerror(status);
    }

    if ((map=(float*)malloc(naxes[0]*naxes[1]*sizeof(float)))==NULL) {
      sprintf(err,"can't allocate image!?");
      complain(err);
    }

    if (newfile) {
      index=naxes[0]*naxes[1];
      // zero
      while (index--) map[index]=0.0;
    } else {
      fits_read_img (ff,TFLOAT,1,naxes[0]*naxes[1],NULL,map,NULL,&status);
      if (status) printerror(status);
    }

    {
      ray aray;
      int x,y;
      float pixsize=side/dim;
      float offset=side/2.0;
      int   idumi=-1;
      float weight;

      if (add_sky_background || add_readout_noise || (poisson_weight!=1)) {
	long iduml1=-1;
	long n;
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv,&tz);
	iduml1=tv.tv_sec+tv.tv_usec;
	n = (iduml1%4096);
	while (n--) poidev(0.5,&idumi);
      }

      // open input

      while (fread(&aray,sizeof(ray),1,stdin)) {
	if (invert) {
	  x=(-aray.p.x+offset)/pixsize;
	  y=(-aray.p.y+offset)/pixsize;
	} else {
	  x=(aray.p.x+offset)/pixsize;
	  y=(aray.p.y+offset)/pixsize;
	}

	if (x<0) continue;
	if (y<0) continue;
	if (x>=dim) continue;
	if (y>=dim) continue;

	if (poisson_weight!=1) {
	  weight=poidev(poisson_weight,&idumi);
	} else {
	  weight=1;
	}

	if (invert) {
	  if (kscale) {
	    map[y*dim+x] -= (modulus(&aray.k)*weight);
	  } else {
	    map[y*dim+x] -= weight;
	  }
	} else {
	  if (kscale) {
	    map[y*dim+x] += (modulus(&aray.k)*weight);
	  } else {
	    map[y*dim+x] += weight;
	  }
	}
      }
      // input is complete
      if (add_sky_background || add_readout_noise) {
	long addr=naxes[0]*naxes[1];
	
	fprintf(stderr,"will add %f electrons sky (poisson) + %f electrons (gaussian) noise\n",poisson_mu,readnoise_per_pixel);
	while (addr--) {
	  if (invert) {
	    if (add_sky_background==1) 
	      map[addr]-=poidev(poisson_mu,&idumi);
	    if (add_readout_noise==1)  
	      map[addr]-=readnoise_per_pixel*gasdev(&idumi);
	  } else {
	    if (add_sky_background==1) 
	      map[addr]+=poidev(poisson_mu,&idumi);
	    if (add_readout_noise==1)  
	      map[addr]+=readnoise_per_pixel*gasdev(&idumi);
	  }
	}
      }
    }

    if (sqroot) {
      int h;
      h=naxes[0]*naxes[1];
      while (h--) map[h]=sqrt(map[h]);
    }

    if (fits_write_img(ff,TFLOAT,1,naxes[0]*naxes[1],map,&status))
      printerror(status);
    if (fits_close_file(ff,&status))
      printerror(status);
  }
  exit(0);
}

void printerror( int status)
{
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
 
  char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
 
  if (status)
    fprintf(stderr, "\n*** Error occurred during program execution ***\n");
 
  fits_get_errstatus(status, status_str);   /* get the error description */
  fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);
 
  /* get first message; null if stack is empty */
  if ( fits_read_errmsg(errmsg) )
    {
      fprintf(stderr, "\nError message stack:\n");
      fprintf(stderr, " %s\n", errmsg);
 
      while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
	fprintf(stderr, " %s\n", errmsg);
    }
 
  exit( status );       /* terminate the program, returning error status */
}      

