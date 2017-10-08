#include "ray.h"
#include "f_2d_tree.h"
#include "zernike.h"

#include "multilayer.h"

// consider structural changes in this software so that this
// #define and #include can be removed...
#define BI_CCD_FUNC
#include "bi_ccd.h"

#include "numrec/include/nr.h"
#include "numrec/include/nrutil.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <argz.h>
#include <time.h>
#include <sys/time.h>

// #define ML_SPEEDUP

float refractivity_mod (float pressure,float temperature);
float beta_t           (float temperature);

void  asp_grad (vec *v,vec *n);
double asp_surf (vec *v);
float n_air_cauchy (float lamda);
float n_air (float lamda);
float n_SiO2_HO2 (float lamda); // fit to 5 values provided in document-2191, allegedly from handbook of optics, vol. 2. hidden option for now
float n_SiO2 (float lamda);
float n_CaF2 (float lamda);
float n_N_BALF4 (float lamda);
float n_N_LLF6 (float lamda);
float n_Ge (float lamda); // dummy placeholder returning n=4 for Ge (IR)
// float (*n_air)(float lamda)=n_air_new;

float (*refractive_index[])(float lamda)={n_SiO2,n_CaF2,n_N_LLF6,n_N_BALF4,n_SiO2_HO2,n_air,n_Ge};


// later these enums should be merged with a common package that is shared
// with multilayer.o library..
typedef enum { asp_SiO2,
	       asp_CaF2,
	       asp_LLF6,
	       asp_N_BALF4,
	       asp_SiO2_HO2,
	       asp_AIR,
	       asp_Ge } MAT;

char *mat_name[]={"SiO2",
		  "CaF2",
		  "LLF6",
		  "N_BALF4",
		  "SiO2_HO2",
		  "air",
		  "Ge" };

typedef struct asph_runpars {
  int flip_sense;
  int refractive_interface;
  int entering;
  int inside;
  int imprint_optic_coords;
  float axial_direction;
  MAT material;
} asph_runpars;

typedef struct airpars {
  float pressure;
  float temperature;
  float pp_H2O;
} airpars;

#ifdef ML_SPEEDUP
typedef struct ml_speedup {
  int   nwave,ntheta;
  float wavemin,wavemax,thetamin,thetamax;
  float *T_ave_LUT;
  float *R_ave_LUT;
} ml_speedup;
#endif

typedef struct multilayerspec {
  char *layerspec;
  optcon *oc;
  float  *oc_const_n;
  float  *layer_thickness;
  multilayer *ml;
  multilayer *ml_reverse;
#ifdef ML_SPEEDUP
  multilayer_speedup *mls;
  multilayer_speedup *mls_reverse;
#endif
} multilayerspec;

#define MAX_ASPHERIC 30

typedef struct asphere {
  double C;
  double k;
  airpars air_pars;
  double A[MAX_ASPHERIC];
  float  r_i;
  float  r_o;
  float  a;
  treelink_control *tlc;
  int    nzernike;
  float  *zernike_amplitudes;
  multilayerspec mlspec;
} asp;


long iduml=-1L;


asp *as;

int surface_interaction(ray *aray,int reflect,asph_runpars *ar,int side);
extern char *usage_str="asphere -R <curv.rad[mm]> -k <conic_const> -a[1-9] <A[A-J] coeff> -ro <r_outer> -ri <r_inner> [-P (print only)}\n";

int main(int argc,char *argv[]) {
  asph_runpars asph_rpars;
  int max_transforms=100;
  int n_transform=0;
  char errstr[512],transform_strings[max_transforms][1024];
  int ind,print_only=0,lens=0,print_dotprod=0,print_interaction_point=0;
  int interference_filter=0;
  lin_trans tf[max_transforms];
  float if_lim[2];
  double r_outer=4300,r_inner=1700;
  asp asph[2];
  asp *asph_p[2];
  char zlist[2048];
  int do_radial_variation=0;
  int do_azimuthal_variation=0;
  float radvar_pv,radvar_mean,azvar_pv,azvar_mean,radvar_Ro,azvar_Ro,var_prob;


  {
    // initialize the run parameters structure
    asph_rpars.flip_sense           = 0;
    asph_rpars.refractive_interface = 0;
    asph_rpars.entering             = 0;
    asph_rpars.imprint_optic_coords = 0;
    asph_rpars.axial_direction      = 0.0;
    asph_rpars.material=asp_SiO2;
  }

  // on the definition of the zernike coefficients
  z_edge_amplitudes=1;

  {
    int i;
    // initialize both structures in this loop
    for (i=0;i<2;i++) {
      asph_p[i]=&asph[i];
      asph[i].air_pars.pressure=760.0;
      asph[i].air_pars.temperature=15.0;
      asph[i].air_pars.pp_H2O=8.0;
      asph[i].nzernike=0;
      asph[i].C=1/1000.0;
      asph[i].k=0;
      asph[i].r_i = r_inner;
      asph[i].r_o = r_outer;
      // zero out treelink business
      asph[i].tlc=NULL;
      // zero out multilayer business
      asph[i].mlspec.layerspec=NULL;
      asph[i].mlspec.ml=NULL;
      asph[i].mlspec.ml_reverse=NULL;
      // zero out aspheric terms
      for (ind=0;ind<MAX_ASPHERIC;ind++) {
	asph[i].A[ind]=0.0;
      }
    }
  }


  as=&asph[0];

  while (--argc) {
    argv++;
    switch (argv[0][0]) {
    case '-':
      switch(argv[0][1]) {
      case 'f':
	asph_rpars.flip_sense=1; // this will change the ordering of interfaces 1 & 2 (if lens) - does nothing if it's not a lens.
	break;
      case 'i':
	print_interaction_point=1;
	break;
      case 'n':
	print_dotprod=1;
	break;
      case 'S':
	asph_rpars.refractive_interface=1;
	switch (argv[0][2]) {
	case '1':
	  // exit; refractive -> (vacuum|air)
	  asph_rpars.entering=0;
	  break;
	case '0':
	default:
	  asph_rpars.entering=1;
	  // entrance; (vacuum|air) -> refractive
	  break;
	}
	break;
      case 'L':
	lens=1;
	break;
      case 'P':
	print_only=1;
	break;
      case 'I':
	interference_filter=1;
	--argc;++argv;	if_lim[0]=atof(argv[0]);
	--argc;++argv;	if_lim[1]=atof(argv[0]);
	break;
      case 'v':
	switch (argv[0][2]) {
	case 'r':
	  // impart a PV variation in 'transmittal' 
	  // from ri to ro on this optic. specify 2 numbers: mean
	  // reflectivity (0-1) and a PV reflectivity (0-1)
	  // since this is a kludge we're not checking for argument count etc.
	  do_radial_variation=1;
	  --argc;++argv;	  radvar_mean=atof(argv[0]);
	  --argc;++argv;	  radvar_pv=atof(argv[0]);
	  break;
	case 'z':
	  // impart a PV variation in 'transmittal' 
	  // that varies with azimuth on this optic. specify 2 numbers: mean
	  // reflectivity (0-1) and a PV reflectivity (0-1)
	  // since this is a kludge we're not checking for argument count etc.
	  // variation should follow: 
	  // reflectivity = azvar_mean + (azvar_pv/2) * cos(az);
	  do_azimuthal_variation=1;
	  --argc;++argv;	  azvar_mean=atof(argv[0]);
	  --argc;++argv;	  azvar_pv=atof(argv[0]);
	  break;
	default:
	  // do nothing
	  continue;
	}
	break;
      case 'l':
	asph_rpars.imprint_optic_coords=1;
	break;
      case 'M':
	--argc; ++argv;
	if (!strcmp(argv[0],"SiO2")) {
	  asph_rpars.material=asp_SiO2;
	} else if (!strcmp(argv[0],"CaF2")) {
	  asph_rpars.material=asp_CaF2;
	} else if (!strcmp(argv[0],"LLF6")) {
	  asph_rpars.material=asp_LLF6;
	} else if (!strcmp(argv[0],"N_BALF4")) {
	  asph_rpars.material=asp_N_BALF4;
	} else if (!strcmp(argv[0],"AIR")) {
	  asph_rpars.material=asp_AIR;
	} else if (!strcmp(argv[0],"Ge")) {
	  asph_rpars.material=asp_Ge;
	} else {
	  // unknown material.
	  fprintf(stderr,"unknown material specification. "
		  "acceptable values are:\n"
		  "\tSiO2\n"
		  "\tCaF2\n"
		  "\tLLF6\n"
		  "\tN-BALF4\n"
		  "\tAIR\n"
		  "\tGe\n"
		  "\nexiting..\n\n");
	  exit(1);
	}
	break;
      case 'T':
	if (n_transform==max_transforms) 
	  complain("too many transforms. increase max_transforms in asphere.c");
	parse_transform_subargs(transform_strings[n_transform++],&argc,&argv);
	break;
      case 'r':
	switch(argv[0][2]) {
	case 'o':
	  --argc; ++argv;
	  asph_p[0]->r_o = asph_p[1]->r_o = r_outer = atof(argv[0]);
	  break;
	case 'i':
	  --argc; ++argv;
	  asph_p[0]->r_i = asph_p[1]->r_i = r_inner = atof(argv[0]);
	  break;
	default:
	  complain("");
	  break;
	}
	break;
      case 'R':      case 'k':      case 'A':      case 'z':      case 'Z':
      case 'p':      case 't':      case 'h':      case 'm':      case 'd':
	switch(argv[0][1]) {
	case 'A':
	  break;
	default:
	  switch(argv[0][2]-'0') {
	  case 1:
	  case 0:
	    as=&asph[0];
	    break;
	  case 2:
	    as=&asph[1];
	    break;
	  default:
	    as=&asph[0];
	    break;
	  }
	  break;
	}
        if (argc<2) {
	  sprintf(errstr,"arguments expected after -%c\n",argv[0][1]);
	  complain(errstr);
	}
	switch(argv[0][1]) {
	case 'm':
	  --argc; ++argv;
	  as->mlspec.layerspec = argv[0];
	  {
	    // allocate multilayer structures for light passing in each 
	    // .. this is the easiest way to operate without redefining
	    // structures & classes.
	    multilayer *mp;
	    mp=(multilayer*)malloc(sizeof(multilayer));
	    mp->ip = NULL;	    mp->mp = NULL;	    mp->lp = NULL;
	    as->mlspec.ml         = mp;
	    mp=(multilayer*)malloc(sizeof(multilayer));
	    mp->ip = NULL;	    mp->mp = NULL;	    mp->lp = NULL;
	    as->mlspec.ml_reverse = mp;
	    mp=NULL;
	  }
	  // parse and initialize the multilayer
	  // using as->mlspec.layerspec
	  {
	    int  nlayer;
	    char *argz;
	    size_t argz_len;
	    char *entry=NULL;
	    char *zlist=as->mlspec.layerspec;
	    int  tmp_argc;

	    argz_create_sep(zlist,' ',&argz,&argz_len);
	    // first count up the number of layers
	    nlayer=0;
	    while ((entry=argz_next(argz,argz_len,entry))) 
	      if ((entry[0]=='-') && (entry[1]=='l')) nlayer++;
	    // allocate.
	    as->mlspec.oc=(optcon*)malloc((nlayer+2)*sizeof(optcon));
	    as->mlspec.oc_const_n=(float*)malloc((nlayer+2)*sizeof(float));
	    as->mlspec.layer_thickness=(float*)malloc(nlayer*sizeof(float));
	    nlayer=0;
	    tmp_argc=argz_count(argz,argz_len);
	    while ((entry=argz_next(argz,argz_len,entry))) {
	      switch(entry[0]) {
	      case 0:
		break;
	      case '-':
		switch(entry[1]) {
		case 'l':
		  if (tmp_argc<=2) 
		    complain("expecting 2 arguments to follow -l");
		  tmp_argc--;
		  entry=argz_next(argz,argz_len,entry);tmp_argc--;
		  as->mlspec.oc[nlayer+1]=get_material(entry);
		  if (as->mlspec.oc[nlayer+1] == const_n_func) {
		    if (tmp_argc <= 2) 
		      complain("expecting 2 arguments to follow -l const");
		    entry=argz_next(argz,argz_len,entry);tmp_argc--;
		    as->mlspec.oc_const_n[nlayer+1]=atof(entry);
		  }
		  entry=argz_next(argz,argz_len,entry);tmp_argc--;
		  as->mlspec.layer_thickness[nlayer]=atof(entry);
		  nlayer++;
		  break;
		default:
		  fprintf(stderr,"huh? %s\n",entry);
		  complain("unexpected switch in multilayer specification.");
		  break;
		}
		break;
	      default:
		fprintf(stderr,"2 -- huh? %s\n",entry);
		complain("unexpected entry in multilayer specification.");
		break;
	      }
	    }
	    // at this point oc, oc_const_n and layer_thickness should be 
	    // populated. now populate the multilayer object
	    // there's a problem with the heterogeneity for initial and final
	    // mediums - routines here are not as flexible as those used in 
	    // the raytrace. the following connections are used only for 
	    // multilayer calculations (reflectivities & transmissions..
	    // fix this at some time so that "air" means the final medium 
	    // and that "SiO2" means the substrate material... etc.
	    if (as == &asph[0]) {
	      // as == &asph[0]
	      as->mlspec.oc[0]                = get_material("air");
	      as->mlspec.oc[nlayer+1]         = get_material("SiO2");
	      as->mlspec.oc_const_n[0]        = 1.0;
	      as->mlspec.oc_const_n[nlayer+1] = 1.5;
	    } else {
	      // as == &asph[1]
	      as->mlspec.oc[0]                = get_material("SiO2");
	      as->mlspec.oc[nlayer+1]         = get_material("air");
	      as->mlspec.oc_const_n[0]        = 1.5;
	      as->mlspec.oc_const_n[nlayer+1] = 1.0;
	    }
	    init_multilayer(as->mlspec.ml,
			    nlayer,
			    as->mlspec.layer_thickness,
			    as->mlspec.oc,
			    as->mlspec.oc_const_n);
	    {
	      // reverse the ordering and the materials list to initialize
	      // the reverse multilayer. reorder:
	      // layer_thickness[nlayer]
	      // oc[nlayer+2]
	      // oc_const_n[nlayer+2]
	      // do this without altering contents of as->mlspec
	      int tmp_ix;
	      optcon* rev_oc=(optcon*)malloc((nlayer+2)*sizeof(optcon));
	      float*  rev_oc_const_n=(float*)malloc((nlayer+2)*sizeof(float));
	      float*  rev_layer_thickness=(float*)malloc(nlayer*sizeof(float));
	      for (tmp_ix=0;tmp_ix<nlayer+2;tmp_ix++) {
		rev_oc[tmp_ix] = as->mlspec.oc[nlayer+2-tmp_ix];
		rev_oc_const_n[tmp_ix] = as->mlspec.oc_const_n[nlayer+2-tmp_ix];
	      }
	      for (tmp_ix=0;tmp_ix<nlayer;tmp_ix++) {
		rev_layer_thickness[tmp_ix]=
		  as->mlspec.layer_thickness[nlayer-tmp_ix];
	      }
	      init_multilayer(as->mlspec.ml_reverse,
			      nlayer,
			      as->mlspec.layer_thickness,
			      as->mlspec.oc,
			      as->mlspec.oc_const_n);
	      free(rev_oc);
	      free(rev_oc_const_n);
	      free(rev_layer_thickness);
	    }
	  }
	  break;
	case 'p':
	  --argc; ++argv;
	  as->air_pars.pressure=760*atof(argv[0]); //convert to torr
	  break;
	case 't':
	  --argc; ++argv;
	  as->air_pars.temperature=atof(argv[0]); // celsius
	  break;
	case 'h':
	  --argc; ++argv;
	  as->air_pars.pp_H2O=atof(argv[0]); // torr
	  break;
	case 'd':
	  --argc; ++argv;
	  // this switch specifies a f_2d_tree type distortion map
	  // that can be interpolated efficiently via a tree algorithm.
	  // input is only a file name for now
	  {
	    if ((as->tlc=(treelink_control*)malloc(sizeof(treelink_control)))
		==NULL) {
	      fprintf(stderr,"can't allocate treelink_control. exiting..\n");
	      exit(1);
	    }
	    as->tlc->treelink_file=argv[0];
	    as->tlc->inited=0;
	    as->tlc->lvs=4;
	    as->tlc->tl=NULL;
	  }
	  break;
	case 'Z':
	  --argc; ++argv;
	  // read in colon delimited list of zernike amplitudes
	  strncpy(zlist,argv[0],2048);
	  {
	    char *argz;
	    size_t argz_len;
	    int i;
	    char *entry=NULL;
	    int nelem;

	    argz_create_sep(zlist,':',&argz,&argz_len);
	    nelem=argz_count(argz,argz_len);
	    as->nzernike=nelem;
	    if ((as->zernike_amplitudes=(float*)malloc(nelem*sizeof(float)))==NULL)
	      complain("cant allocate zernike_amplitudes..\n");
	    for (i=0;i<nelem;i++) {
	      entry=argz_next(argz,argz_len,entry);
	      as->zernike_amplitudes[i]=atof(entry);
	    }
	  }
	  break;
	case 'z':
	  --argc; ++argv; 
	  (as->A)[0]=atof(argv[0]);
	  break;
	case 'R':
	  --argc; ++argv; 
	  as->C=1.0/atof(argv[0]);
	  break;
	case 'k':
	  --argc; ++argv; 
	  as->k=atof(argv[0]);
	  break;
	case 'A':
	  ind=atoi(argv[0]+2);
	  if (ind<1 || ind>MAX_ASPHERIC-1) {
	    complain("huh..?");
	  }
	  --argc; ++argv; 
	  (as->A)[ind]=atof(argv[0]);
	  
	  (as->A)[ind] *= 1e3/pow(1e3,ind); // using Lynn Sepalla's units
	  // now (as->A)[ind] has units of mm/mm^{ind} - Lynn's units were
	  // m/mm^{ind}.
	  break;
	default:
	  break;
	}
        break;
      default:
        sprintf(errstr,"unknown switch: %s\n",argv[0]);
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
    // initialize the random number generator used here
    struct timeval tv;
    struct timezone tz;
    int i;
    gettimeofday(&tv,&tz);
    iduml = tv.tv_sec + tv.tv_usec;
    i=iduml%4096;
    while (--i) ran2(&iduml);
    // iduml is initialized with a seed originating from the time
  }

  if (do_radial_variation) {
    // compute these numbers here.
    radvar_Ro = radvar_mean / 
      ( 1 + radvar_pv/(r_outer-r_inner)*
	(r_inner - (2./3.)*((pow(r_outer,3)-pow(r_inner,3))
			    /(pow(r_outer,2)-pow(r_inner,2)))));
    radvar_pv /= radvar_Ro;
    // ready to use radvar_Ro and radvar_pv
    fprintf(stderr,"radvar pars: Ro = %f pv=%f\n",radvar_Ro,radvar_pv);
  }
  if (do_azimuthal_variation) {
    azvar_Ro=azvar_mean;
    // ready to use azvar_Ro and azvar_pv
    fprintf(stderr,"azvar pars: Ro = %f pv=%f\n",azvar_Ro,azvar_pv);
  }

  as=&asph[0];
  as->a=r_outer;
  {
    char *s;

    if (lens && asph_rpars.refractive_interface) {
      complain("shouldn't specify both <single> refractive interface and <lens> at the same time. parsing ambiguity will result!!\n");
    }

    if (lens) {
      fprintf(stderr,
	      "             Mat.  = %s\n",mat_name[asph_rpars.material]);
      s="(1)";
    } else {
      s="   ";
    }

    if (asph_rpars.refractive_interface) {
      fprintf(stderr,
	      "             Mat.  = %s\n",mat_name[asph_rpars.material]);
      if (asph_rpars.entering) {
	fprintf(stderr,
		"             (entering)\n");
      } else {
	fprintf(stderr,
		"             (exiting)\n");
      }
      s="   ";
    }
  
    fprintf(stderr,
	    "asphere setup:  \n"
	    "%s           R     = %+g mm\n"
	    "%s           k     = %+g\n",
	    s,1.0/as->C,
	    s,as->k);

    if (lens || asph_rpars.refractive_interface) {
      fprintf(stderr,
	      "%s           p_atm[torr] = %+g\n",s,as->air_pars.pressure);
      fprintf(stderr,
	      "%s           T[C]        = %+g\n",s,as->air_pars.temperature);
      fprintf(stderr,
	      "%s           pp_H2O[torr]= %+g\n",s,as->air_pars.pp_H2O);
    }
    if (as->mlspec.layerspec) {
      int i;
      if (1)
	fprintf(stderr,
		"%s           multilayer  = %s\n",s,as->mlspec.layerspec);
      fprintf(stderr,"%s              got a total of %d layers..\n",s,as->mlspec.ml->nlayer);
      if (1) {
	fprintf(stderr,"%s              initial medium %p\n",s,as->mlspec.oc[0]);
	for (i=0;i<as->mlspec.ml->nlayer;i++) {
	  fprintf(stderr,"%s              oc %p thickness %f\n",
		  s,as->mlspec.oc[i+1],as->mlspec.layer_thickness[i]);
	}
	fprintf(stderr,"%s              final medium %p\n",s,as->mlspec.oc[as->mlspec.ml->nlayer+1]);
      }

    }
    for (ind=0;ind<MAX_ASPHERIC;ind++) {
      if ((as->A)[ind]!=0) 
	fprintf(stderr,"%s           A[%02d] = %+g mm^{%d}\n",s,ind,(as->A)[ind],1-ind);
    }

    fprintf(stderr,"%s           tlc   = %p\n",s,as->tlc);

    if (as->tlc!=NULL) {
      fprintf(stderr,"%s              file = %s\n",s,as->tlc->treelink_file);
    }

    if (as->nzernike) {
      int j,n,m;
      if (z_edge_amplitudes == 1) {
	fprintf(stderr,"%s zernike amplitudes at aperture edge [not coefficients]:\n",s);
      } else {
	fprintf(stderr,"%s zernike coefficients [not amplitudes]:\n",s);
      }
      for (j=0;j<as->nzernike;j++) {
	n=ceil((-3+sqrt(9+8*j))/2);
	m=2*j-n*(n+2);
	fprintf(stderr,"%s (j =%3d n =%3d m =%3d) : %8.4g [microns]\n",
		s,j,n,m,as->zernike_amplitudes[j]*1e3);
      }
    }
    
    if (lens) {
      as=&asph[1];
      as->a=r_outer;

      s="(2)";
      fprintf(stderr,
	      "%s           R     = %+g mm\n"
	      "%s           k     = %+g\n",
	      s,1.0/as->C,s,as->k);
      fprintf(stderr,
	      "%s           p_atm[torr] = %+g\n",s,as->air_pars.pressure);
      fprintf(stderr,
	      "%s           T[C]        = %+g\n",s,as->air_pars.temperature);
      fprintf(stderr,
	      "%s           pp_H2O[torr]= %+g\n",s,as->air_pars.pp_H2O);
      if (as->mlspec.layerspec) {
	int i;
	if (1)
	  fprintf(stderr,
		  "%s           multilayer  = %s\n",s,as->mlspec.layerspec);
	fprintf(stderr,"%s              got a total of %d layers..\n",s,as->mlspec.ml->nlayer);
	if (1) {
	  fprintf(stderr,"%s              initial medium %p\n",s,as->mlspec.oc[0]);
	  for (i=0;i<as->mlspec.ml->nlayer;i++) {
	    fprintf(stderr,"%s              oc %p thickness %f\n",
		    s,as->mlspec.oc[i+1],as->mlspec.layer_thickness[i]);
	  }
	  fprintf(stderr,"%s              final medium %p\n",s,as->mlspec.oc[as->mlspec.ml->nlayer+1]);
	} 
      }

      for (ind=0;ind<MAX_ASPHERIC;ind++) {
	if ((as->A)[ind]!=0) 
	  fprintf(stderr,"%s           A[%02d] = %+g mm^{%d}\n",s,ind,(as->A)[ind],1-ind);
      }

      fprintf(stderr,"%s           tlc   = %p\n",s,as->tlc);

      if (as->tlc!=NULL) {
	fprintf(stderr,"%s              file = %s\n",s,as->tlc->treelink_file);
      }

      if (as->nzernike) {
	int j,n,m;
	fprintf(stderr,"%s zernike amplitudes:\n",s);
	for (j=0;j<as->nzernike;j++) {
	  n=ceil((-3+sqrt(9+8*j))/2);
	  m=2*j-n*(n+2);
	  fprintf(stderr,"%s (j =%3d n =%3d m =%3d) : %8.4g [microns]\n",
		  s,j,n,m,as->zernike_amplitudes[j]*1e3);
	}
      }
    }
  }

  {
    if (asph_rpars.flip_sense && lens) {
      fprintf(stderr,"flipping ordering of surfaces for this lens: 1->2; 2->1\n");
      asph_p[0]=&asph[1];
      asph_p[1]=&asph[0];
    }
    if (asph_rpars.flip_sense && asph_rpars.refractive_interface) {
      asph_rpars.entering = !asph_rpars.entering;
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

  if (print_only) {
    vec v;
    v.z=v.y=0;
    
    as=&asph[0];
    for (v.x=-r_outer;v.x<=r_outer;v.x+=r_outer/1000) {
      if (fabs(v.x)>r_inner)
	printf("%g %g\n",v.x,asp_surf(&v));
    }
    if (lens) {
      as=&asph[1];
      for (v.x=r_outer;v.x>=-r_outer;v.x-=r_outer/1000) {
	if (fabs(v.x)>r_inner)
	  printf("%g %g\n",v.x,asp_surf(&v));
      }
      v.x=-r_outer;
      as=&asph[0];
      printf("%g %g\n",v.x,asp_surf(&v));
    }
    exit(0);
  }
  {
    ray aray;
    vec normal;
    //    vec normal2;
    vec inray_p,outray_p;
    double dp1,dp2;
    float r;
    int echo_output_n=0;
    /* 
    float lamda;
    for (lamda=400;lamda<=1100;lamda+=5) {
      printf("%f %f %f %f %f\n",
	     lamda,
	     (*refractive_index[SiO2])(lamda),
	     (*refractive_index[CaF2])(lamda),
	     (*refractive_index[LLF6])(lamda),
	     (*refractive_index[N_BALF4])(lamda));
    }
    exit(0);
    */

    if (lens) {
      // positive if apex of 2nd surface is more positive than the apex of 1st
      asph_rpars.axial_direction = (asph_p[1]->A[0] - asph_p[0]->A[0]);
    } else {
      // not really used. set to 1 in case.
      asph_rpars.axial_direction = 1.0;
    }

    while (fread(&aray,sizeof(ray),1,stdin)) {
      if (asph_rpars.imprint_optic_coords) echo_output_n=1;
      if (do_azimuthal_variation || do_radial_variation) var_prob = 1;

      if (n_transform>0) {
	int i=0;
	while (i<n_transform) {
	  tran_ray(&aray,&tf[i],(tf[i].invert==0)?1:-1);
	  i++;
	}
      }
      
      // compute dp1 here in case incidence angles are required.
      if (print_dotprod) {
	ray_surface_collapse(&aray,asp_surf);
	asp_grad(&aray.p,&normal);
	dp1=fabs(dot_prod(&aray.k,&normal)/modulus(&aray.k));
      }
      
      // here work out the loop in a modular way
      
      if (0) { // -vz and -vr switches are disabled here..
	int nasp,ai,reflect,reflect_result;

	if (asph_rpars.refractive_interface || lens) {
	  reflect=0;
	  if (lens) 
	    nasp=2;
	  else
	    nasp=1;
	} else {
	  // do a simple reflection
	  reflect=1;
	  nasp=1;
	}

	ai=0;
	asph_rpars.inside=0;

	while (ai<nasp) {
	  as=asph_p[ai];
	  reflect_result = surface_interaction(&aray,reflect,&asph_rpars,ai);
	  switch(reflect_result) {
	  case -1:
	    // destroyed ray.
	    goto skip_ray;
	    break;
	  case 0:
	    // a transmission.
	    if (reflect==0) {
	      if ((ai==0) && (asph_rpars.inside==1)) 
		goto skip_ray;
	      // transmission, just as asked. so far so good..
	      if (asph_rpars.inside)
		asph_rpars.inside=0;
	      else
		asph_rpars.inside=1;
	    } else {
	      // reflect==1 should override possibility of transmission.
	      // skip this.. 
	      goto skip_ray;
	    }
	    break;
	  case 1:
	    // a reflection.
	    // do not alter value of asph_rpars.inside.
	    if (reflect==0) {
	      // reflection while asking for transmission.
	      // internal reflections computed in this case. test ai etc.
	      if (ai==1) {
		// on the 2nd surface.
		// internal reflection. 
		ai -= 2;
	      } else {
		// on the 1st surface.
		if (asph_rpars.inside==1) {
		  // do nothing, ray is still in play.
		} else {
		  // reflection off the first surface from the outside.
		  // ray is dead..
		  goto skip_ray;
		}
	      }
	    } else {
	      // reflection, just as asked. so far so good..
	    }
	    break;
	  default:
	    fprintf(stderr,
		    "surface_interaction unexpected return..\nexiting..\n");
	    exit(1);
	    break;
	  }
	  ai++;
	}
	goto output_ray;
      } else {
	as=&asph[0];
	as=asph_p[0];
	ray_surface_collapse(&aray,asp_surf);

	memcpy(&inray_p,&aray.p,sizeof(vec));

	r=sqrt(pow(aray.p.y,2.0)+pow(aray.p.x,2.0));
	if ((r<r_outer) && (r>r_inner)) {
	  if (do_azimuthal_variation || do_radial_variation) {
	    float theta=atan2(aray.p.y,aray.p.x);    // [-PI..+PI]
	    float rho=(r-r_inner)/(r_outer-r_inner); // [0..1]
	    if (do_azimuthal_variation) {
	      // this change turns azimuthal variation into a lateral..
	      var_prob *= azvar_Ro*(1 + 0.5*azvar_pv*(r/r_outer)*cos(theta) );
	    }
	    if (do_radial_variation) {
	      var_prob *= radvar_Ro*(1 - radvar_pv*rho);
	    }
	    //	  fprintf(stderr,"comparing var_prob = %g\n",var_prob);
	    if (rand() > RAND_MAX*var_prob) 
	      goto skip_ray;
	  }
	  if (asph_rpars.imprint_optic_coords) {
	    // compute optic coords for this ray.
	    float theta=atan2(aray.p.y,aray.p.x);    // [-PI..+PI]
	    float rho=(r-r_inner)/(r_outer-r_inner); // [0..1]
	    if ((0.5-fabs(0.5
			  -(theta/(2*M_PI/24)-floor(theta/(2*M_PI/24))))<0.025)
		||
		(0.5-fabs(0.5
			  -(rho/(1.0/5)-floor(rho/(1.0/5))))<0.025)) {
	      echo_output_n*=2;
	    }
	  }
	  // ray is reflected. assume unity reflectivity for now.
	  // surface_normal(&aray.p,asp_surf,&normal);
	  asp_grad(&aray.p,&normal);

	  if (interference_filter && !(lens || asph_rpars.refractive_interface)) {

	    fprintf(stderr,"not handling interference filter for a reflector properly..\n\tfix this soon..\n");
	    float a = 2*M_PI/fabs(dot_prod(&aray.k,&normal));
	    if ((a>if_lim[1]) || (a<if_lim[0])) goto skip_ray;
	  }

	  dp1=fabs(dot_prod(&aray.k,&normal)/
		   modulus(&aray.k));

	  if (as->mlspec.layerspec != NULL){
	    // compute the multilayer thing.
	    //	    float theta,lambda,kmod;
	    //	    theta=acos(dp1);
	    //	    kmod=modulus(&aray.k);
	    //	    lamda=0.1*2*M_PI/kmod;
	    //	    compute_multilayer(theta,lambda,as->mlspec.ml);
	    // if this is not a lens can go ahead and test to see
	    // if ray should be skipped.
	    //	    if (ran() > as->mlspec.ml->R_ave) {
	    //	      process=0;
	    //	    }
	    //	  } else {
	    //	    process=1;
	  }


	  if (lens || asph_rpars.refractive_interface) {
	    float lamda;
	    float kmod;

	    kmod=modulus(&aray.k);
	    lamda=0.1*2*M_PI/kmod; // nm

	    // apply snell's law etc.
	    if (1) {
	      if (asph_rpars.refractive_interface) {
		if (asph_rpars.material==asp_AIR) { // vacuum -> air case..

		  if (interference_filter) {
		    fprintf(stderr,"not handling interference filter properly "
			    "for VAC->AIR interface (refractive_interface==1,ma"
			    "terial==AIR)\n\tfix this soon..\n");
		  }

		  {
		    float n_fr,n_to;
		    n_fr = ((asph_rpars.entering)?1:(*refractive_index[asph_rpars.material])(lamda));
		    n_to = ((asph_rpars.entering)?(*refractive_index[asph_rpars.material])(lamda):1);
		    if (refract_ray(&aray.k,&normal,n_fr,n_to)) {
		      //		    aray.n = n_to;
		      // send the photon on its way.
		      goto output_ray;
		    } else {
		      goto skip_ray;
		    }
		  }
		} else {

		  if (interference_filter) {
		    if (asph_rpars.entering) {
		      float sin_theta_prime=
			sqrt(1-dp1*dp1)/(*refractive_index[asph_rpars.material])(lamda);
		      float a=(2*M_PI/kmod)/sqrt(1-sin_theta_prime*sin_theta_prime);
		      if ((a>if_lim[1]) || (a<if_lim[0])) goto skip_ray;
		    } else {
		      fprintf(stderr,"not handling interference filter properly "
			      "for 'exiting' refractive interface\n"
			      "\tfix this soon..\n");
		    }
		  } 
		  {
		    float n_fr,n_to;
		    n_fr=((asph_rpars.entering)?n_air(lamda):
			  (*refractive_index[asph_rpars.material])(lamda));
		    n_to=((asph_rpars.entering)?(*refractive_index[asph_rpars.material])(lamda):
			  n_air(lamda));
		    if (refract_ray(&aray.k,&normal,n_fr,n_to)) {
		      //		    aray.n = n_to;
		      // send the photon on its way.
		      goto output_ray;
		    } else {
		      goto skip_ray;
		    }
		  }
		}
	      }

	      // if control arrives here then this is a lens.

	      if (interference_filter) {
		float sin_theta_prime=
		  sqrt(1-dp1*dp1)/(*refractive_index[asph_rpars.material])(lamda);
		float a=(2*M_PI/kmod)/sqrt(1-sin_theta_prime*sin_theta_prime);
		if ((a>if_lim[1]) || (a<if_lim[0])) goto skip_ray;
	      }

	      {
		float n_fr,n_to;

		n_fr=n_air(lamda);
		n_to=(*refractive_index[asph_rpars.material])(lamda);
		if (refract_ray(&aray.k,
				&normal,
				n_fr,
				n_to)) {
		  as=&asph[1];
		  as=asph_p[1];
		  //		aray.n = n_to;
	      
		  ray_surface_collapse(&aray,asp_surf);
	      
		  memcpy(&outray_p,&aray.p,sizeof(vec));
	      
		  r=sqrt(pow(aray.p.y,2.0)+pow(aray.p.x,2.0));
		  if (r<r_outer && r>r_inner) {
		    if (asph_rpars.imprint_optic_coords) {
		      // compute optic coords for this ray.
		      float theta=atan2(aray.p.y,aray.p.x);    // [-PI..+PI]
		      float rho=(r-r_inner)/(r_outer-r_inner); // [0..1]
		      if ((0.5-fabs(0.5
				    -(theta/(2*M_PI/24)-floor(theta/(2*M_PI/24))))<0.025)
			  ||
			  (0.5-fabs(0.5
				    -(rho/(1.0/5)-floor(rho/(1.0/5))))<0.025)) {
			echo_output_n*=2;
		      }
		    }
		
		    // surface_normal(&aray.p,asp_surf,&normal);
		    asp_grad(&aray.p,&normal);
		    if (1) {
		      float tmp;
		      tmp=n_to;
		      n_to=n_fr;
		      n_fr=tmp;
		    } else {
		      n_fr=(*refractive_index[asph_rpars.material])(lamda);
		      n_to=n_air(lamda);
		    }
		    if (refract_ray(&aray.k,&normal,n_fr,n_to)) {
		      //		    aray.n=n_to;
		      dp2=fabs(dot_prod(&aray.k,&normal)/
			       modulus(&aray.k));
		      r=sqrt(pow(aray.p.y,2.0)+pow(aray.p.x,2.0));
		      if (r<r_outer && r>r_inner) {
			goto output_ray;
		      }
		      goto skip_ray;
		    }
		    goto skip_ray;
		  }
		  goto skip_ray;
		}
	      }
	      goto skip_ray;
	    }
	  } else {
	    // reflect off of a mirror.
	    if (dot_prod(&aray.k,&normal)>0) {
	      scalevec(&normal,-1);
	      //	    normal.x *= -1;
	      //	    normal.y *= -1;
	      //	    normal.z *= -1;
	    }
	    reflect_ray(&aray.k,&normal);
	    goto output_ray;
	  }
	} else {
	  // ray is missed
	  goto skip_ray;
	}
      }

    output_ray:
      if (n_transform>0) {
	int i=n_transform;
	while (i>0) {
	  i--;
	  tran_ray(&aray,&tf[i],(tf[i].invert==0)?-1:1);
	}
      }

      //      if (transform_string) tran_ray(&aray,&tf,-1);
      if (print_dotprod || print_interaction_point) {
	if (print_dotprod) {
	  if (lens) {
	    fprintf(stdout,"%f %f\n",180.0/M_PI*acos(dp1),180.0/M_PI*acos(dp2));
	  } else {
	    fprintf(stdout,"%f\n",180.0/M_PI*acos(dp1));
	  }
	}
	if (print_interaction_point) {
	  if (lens) {
	    fprintf(stdout,"%f %f %f %f %f %f\n",inray_p.x,inray_p.y,inray_p.z,
		    outray_p.x,outray_p.y,outray_p.z);
	  } else {
	    fprintf(stdout,"%f %f %f\n",inray_p.x,inray_p.y,inray_p.z);
	  }
	}
      } else {
	fwrite(&aray,sizeof(ray),1,stdout);
	if (asph_rpars.imprint_optic_coords) {
	  while (--echo_output_n) {
	    fwrite(&aray,sizeof(ray),1,stdout);
	  }
	}
      }
    skip_ray:
      continue;
    }  
  }
  exit(0);
}
  
void asp_grad (vec *v,vec *n) {
  double r=sqrt(pow(v->x,2.0)+pow(v->y,2.0));
  vec rhat,zhat,thetahat;
  double a,rscale;
  int    ind;
  float dz_drho,dz_dtheta,theta_scale;

  theta_scale=0;

  rhat.x=v->x;
  rhat.y=v->y;
  rhat.z=0;
  unitvec(&rhat);
  zhat.x=zhat.y=0;
  zhat.z=1;
  unitvec(&zhat);

  cross_prod(&zhat,&rhat,&thetahat);
  unitvec(&thetahat);

  a=sqrt(1-(1+as->k)*pow(as->C*r,2.0));
  rscale=(as->C*r)/(1+a)*(2 + pow(as->C*r,2)*(1+as->k)/(a*(1+a)));
  for (ind=1;ind<MAX_ASPHERIC;ind++) {
    rscale += ind*(as->A)[ind]*pow(r,ind-1);
  }


  double ts1=0,ts2=0,rs1=0,rs2=0;

  if (as->nzernike) {
    zernike_distortions(r/as->a,atan2(v->y,v->x),
			as->zernike_amplitudes,as->nzernike,
			&dz_drho,&dz_dtheta);
    rscale     += dz_drho/as->a;
    theta_scale += dz_dtheta/r;
    rs1     = dz_drho/as->a;
    ts1     = dz_dtheta/r;
  }

  if (as->tlc != NULL) {
    point p;
    float zval,dzdx,dzdy;
    p.x=v->x;    p.y=v->y;
    float theta=atan2(p.y,p.x);

    float costheta=cos(theta);
    float sintheta=sin(theta);

    tl_interp(as->tlc,&p,&zval,&dzdx,&dzdy);
    rscale      += (costheta*dzdx + sintheta*dzdy);
    theta_scale += (-sintheta*dzdx + costheta*dzdy);
    rs2 = (costheta*dzdx + sintheta*dzdy);
    ts2 = (-sintheta*dzdx + costheta*dzdy);
  }

  //  if (as->tlc && as->nzernike) {
    //    fprintf(stderr,"x,y= %g,%g thetascales = %g,%g rscales = %g,%g\n",
    //	    v->x,v->y,ts1,ts2,rs1,rs2);
  //  }
  
  scalevec(&rhat,rscale);
  scalevec(&thetahat,theta_scale);
  vec_add(&rhat,&thetahat,n);
  vec_diff(n,&zhat,n);
  unitvec(n);
}

double asp_surf (vec *v) {
  double r=sqrt(pow(v->x,2.0)+pow(v->y,2.0));
  double z=as->C*pow(r,2.0)/(1+sqrt(1-(1+as->k)*pow(as->C*r,2.0)));
  float  dz_drho,dz_dtheta;
  int    ind;
  for (ind=0;ind<MAX_ASPHERIC;ind++) {
    if (fabs((as->A)[ind])!=0) {
      z+=(as->A)[ind]*pow(r,ind);
    }
  }
  if (as->nzernike)
    z+=zernike_distortions(r/as->a,atan2(v->y,v->x),
			   as->zernike_amplitudes,as->nzernike,
			   &dz_drho,&dz_dtheta);
  if (as->tlc != NULL) {
    point p;
    float zval,dzdx,dzdy;
    p.x=v->x;
    p.y=v->y;
    tl_interp(as->tlc,&p,&zval,&dzdx,&dzdy);
    // sign of zval is not a sure thing yet
    z += zval;
  }
  return(z-v->z);
}

float n_SiO2 (float lamda) {
  // lamda is in nm.
  return(1.45229948
         +2717.00875*pow(lamda,-2)
         +4.88616329E+7*pow(lamda,-4)
         +lamda*(-1.85896965E-06
                 -2.76150231E-10*lamda));
}

float n_SiO2_HO2 (float lamda) {
  // lamda is in nm.
  return(1.50310719+lamda*(-0.000135856943+lamda*(1.27841076E-07+lamda*(-4.47723698E-11))));
}

float n_CaF2 (float lamda) {
  // lamda is in nm.
  float C1  = 0.001719,   C2 = 1.46253,   C3 = 0.111545;
  float C4  =-0.0559683,  C5 = 0.0347141, C6 =-0.00323392;
  float C7  = 0.0181429,  C8 = 9.8555e-05,C9 = 0.0104573;
  float C10 = 0.000217936,C11= 0.002986;
  float lamda_min=0.1311,lamda_max=9.4291;
  float L_lam=1/(1e-6*lamda*lamda-C1);
  float L_min=1/(lamda_min*lamda_min-C1);
  float L_max=1/(lamda_max*lamda_max-C1);
  float deltaL=(L_min-L_max)/2.0;
  float L_ave=(L_max+L_min)/2.0;
  float mu=(L_lam-L_ave)/deltaL;
  float delta_lam=(lamda_max*lamda_max-lamda_min*lamda_min)/2.0;
  float lamda_ave=(lamda_max*lamda_max+lamda_min*lamda_min)/2.0;
  float nu=(1e-6*lamda*lamda-lamda_ave)/delta_lam;
  
  //  fprintf(stderr,"mu=%g nu=%g\n");
  // return resnik formula
  return(C2 + nu*(C4 +nu*(C6 +nu*(C8 +nu*C10))) +
	      mu*(C3 +mu*(C5 +mu*(C7 +mu*(C9 +mu*C11)))));
  // this was an old cheapo formula.
  //  return(1.43170655
  //         +1942.90571*pow(lamda,-2)
  //         +2.59218104E7*pow(lamda,-4)
  //         +lamda*(-8.84210039E-6
  //                 +2.98957457E-9*lamda));
}

float n_Ge (float lambda) {
  // lambda is in nm.
  // return constant dummy value of n=4.
  //  float C1=1.31004128,C2=0.142038259,C3=0.964929351;
  //  float C4=0.0079659645,C5=0.0330672072,C6=109.19732;
  //  float l2=pow(lambda/1e3,2);
  // return sellmeier formula
  //  return(sqrt(1 + C1*l2/(l2-C4) + C2*l2/(l2-C5) + C3*l2/(l2-C6)));
  return(4.0034);
}

float n_N_BALF4 (float lambda) {
  // lambda is in nm.
  float C1=1.31004128,C2=0.142038259,C3=0.964929351;
  float C4=0.0079659645,C5=0.0330672072,C6=109.19732;
  float l2=pow(lambda/1e3,2);
  // return sellmeier formula
  return(sqrt(1 + C1*l2/(l2-C4) + C2*l2/(l2-C5) + C3*l2/(l2-C6)));
}

float n_N_LLF6 (float lambda) {
  // lambda is in nm.
  float C1=1.22932376,C2=0.0755174062,C3=0.988839837;
  float C4=0.00877776451,C5=0.0480429329,C6=112.136098;
  float l2=pow(lambda/1e3,2);
  // return sellmeier formula
  return(sqrt(1 + C1*l2/(l2-C4) + C2*l2/(l2-C5) + C3*l2/(l2-C6)));
}

float n_silica_old (float lamda) {
  // this was based on an older fit.. only different by ~ 1 part in 1e4 
  // compared to the newer one (above)
  // lamda is in nm.
  return(1.45597219
	 +2464.48853*pow(lamda,-2)
	 +55607448.0*pow(lamda,-4)
	 +lamda*(-8.86795351E-06
		 +6.47343068E-10*lamda));
}

float n_air (float lamda) {
  // lamda is in nm.
  // old return values..
  // this is cauchy's formula.
  //  return(1.000287566
  //	 +1.3412e-18*pow(lamda*1e-9,-2)
  //	 +3.777e-32*pow(lamda*1e-9,-4));
  float refractivity,n;
  // this routine uses the global structure air_pars to fine tune n.
  // air_pars is modified in command line arguments.
  // this estimate is from Edlen, 1953.
  float alpha=1/273.15;
  refractivity=64.328
    +29498.1/(146-pow(lamda*1e-3,-2))
    +255.4/(41-pow(lamda*1e-3,-2));
  // scale refractivity according to pressure & temperature estimates.
  // (Barrell 1951)
  refractivity *= refractivity_mod(as->air_pars.pressure,as->air_pars.temperature);
  // and modify according to humidity. (Barrell, 1939)
  refractivity -= ((6.24e-2 - 6.80e-4*pow(lamda*1e-3,-2))*as->air_pars.pp_H2O/
		   (1+alpha*as->air_pars.temperature));
  n=1+1e-6*refractivity;
  return(n);
}

float n_air_cauchy (float lamda) {
  // lamda is in nm.
  // old return values..
  // this is cauchy's formula.
  return(1.000287566
  	 +1.3412e-18*pow(lamda*1e-9,-2)
	 +3.777e-32*pow(lamda*1e-9,-4));
}

float refractivity_mod (float pressure,float temperature) {
  float alpha=1/273.15;
  return((pressure/760)*
	 (1+beta_t(temperature)*pressure)/(1+beta_t(15)*760)*
	 (1+alpha*15)/(1+alpha*temperature));
}

float beta_t (float temperature) {
  return(1.049-0.0157*temperature);
}

int surface_interaction(ray *t_ray,int reflect,asph_runpars *ar,int side) {
  // routine to return result of the interaction between ray and surface.
  // if reflect=1 then reflect off of the surface.
  // otherwise, depending on the contents of memory pointed to by asp *as 
  // the new ray is returned in place, with the return value signifying
  // whether the surface interaction produced a reflection (1) or not (0)
  vec   normal;
  float r,theta,lambda;
  double dp,rn;
  float n_fr,n_to;
  multilayer *mlp;

  ray_surface_collapse(t_ray,asp_surf);
  r=sqrt(pow(t_ray->p.y,2.0)+pow(t_ray->p.x,2.0));
  if ((as->r_o - r)*(r - as->r_i) < 0) {
    fprintf(stderr,"falls outside..\n");
    return(-1); // falls outside annulus
  }
  asp_grad(&t_ray->p,&normal);
  dp=dot_prod(&t_ray->k,&normal)/modulus(&t_ray->k);
  // depending on the sign of the dot product, fix the normal vector
  if (dp>0) {
    scalevec(&normal,-1);
  }

  if (as->mlspec.layerspec != NULL) {
    theta=acos(fabs(dp));
    lambda=0.1*2*M_PI/modulus(&t_ray->k);

    if ((side==0) && ar->inside) {
      mlp=as->mlspec.ml_reverse;
    } else {
      mlp=as->mlspec.ml;
    }

    compute_multilayer(theta,lambda,mlp);
    // roll the die to see whether this will be a reflection or transmission
    // or neither
    rn=ran2(&iduml);

    rn -= mlp->R_ave;
    if (rn < 0) {
      // it's a reflection
      goto refl;
    } else {

      rn -= mlp->T_ave;
      if (rn < 0) {
	// it's a transmission
	goto tran;
      } else {
	// it's an absorption. destroy..
	return(-1);
      }
    }
  }

  if (reflect==1) {
    // do a reflection
    goto refl;
  } else {
    // do a transmission
    goto tran;
  }
  
 refl:
  reflect_ray(&t_ray->k,&normal);
  return(1);
 tran:
  if (as->mlspec.layerspec == NULL) {
    // lambda hasn't been calculated yet. calculate here..
    lambda=0.1*2*M_PI/modulus(&t_ray->k);
  }
  if (ar->refractive_interface) {
    if (ar->material==asp_AIR) { 
      // vacuum -> air interface case..
      n_fr=((ar->entering)?
	    1:
	    (*refractive_index[ar->material])(lambda));
      n_to=((ar->entering)?
	    (*refractive_index[ar->material])(lambda):
	    1);
    } else {
      n_fr=((ar->entering)?
	    n_air(lambda):
	    (*refractive_index[ar->material])(lambda));
      n_to=((ar->entering)?
	    (*refractive_index[ar->material])(lambda):
	    n_air(lambda));
    }
  } else {
    // a lens.
    // pressure, humidity etc. are already set in parameters
    // pointed to by asph_p (used in n_air)..
    if (side == 0) {
      n_fr=n_air(lambda);
      n_to=(*refractive_index[ar->material])(lambda);
    } else {
      n_fr=(*refractive_index[ar->material])(lambda);
      n_to=n_air(lambda);
    }
  }
  if (refract_ray(&t_ray->k,&normal,n_fr,n_to)) {
    return(0);
  } else {
    return(-1);
  }
}
