#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "f_2d_tree.h"

void tl_interp(treelink_control *tlc,point *p,float *zval,float *dzdx,float *dzdy) {
  int ntl=5;
  treelink *tl_list[ntl];
  if (tlc->inited==0)
    init_tlc(tlc);
  int nn=tree_nearest_links(tlc->tl,p,tlc->lvs,tl_list,ntl);
  if (nn>2) { 
    /* should normally be the case, since tree was pruned >2 above */
    /* can do vector calculations maybe. */
    /* with model z(x,y) = mx*x + my*y + b */
    
    /* /x2 xy x\ /mx\   / xz \ */
    /* |xy y2 y| |my| = | yz | */
    /* \x  y  I/ \b /   \ z  / */
    /* inverse of the sampling matrix is the transpose of cofactor matrix */
    /* divided through by the determinant of the sampling matrix. */
    double XX,XY,YY,X,Y,I,XZ,YZ,Z;
    double x,y,z,x0,y0,z0;
    x0=tl_list[0]->pt->x;
    y0=tl_list[0]->pt->y;
    z0=*(tl_list[0]->val);
    int k;
    XX=0.0;    XY=0.0;    YY=0.0;
    X=0.0;     Y=0.0;     I=0.0;
    XZ=0.0;    YZ=0.0;    Z=0.0;
    for (k=0;k<nn;k++) {
      x = tl_list[k]->pt->x - x0;
      y = tl_list[k]->pt->y - y0;
      z = *(tl_list[k]->val)- z0;
      XX += x*x;      XY += x*y;      YY += y*y;
      X  += x;        Y  += y;        I  += 1;
      XZ += x*z;      YZ += y*z;      Z  += z;
      //      fprintf(stderr,"sample %d: x,y,z = %f %f %f\n",(int)I,x+x0,y+y0,z+z0);
    }
    double DET=XX*(YY*I-Y*Y)+XY*(Y*X-I*XY)+X*(XY*Y-X*YY);

    if (DET==0) {
      fprintf(stderr,"youch! DET is zero (x,y)=(%g,%g)..\n",
      	      p->x,p->y);
      goto output_nearest;
    }
    double mx,my,b,mxDET,myDET,bDET;

    mxDET=(YY*I-Y*Y)*XZ+(X*Y-I*XY)*YZ+(XY*Y-X*YY)*Z;
    myDET=(Y*X-I*XY)*XZ+(I*XX-X*X)*YZ+(X*XY-XX*Y)*Z;
    bDET =(XY*Y-YY*X)*XZ+(X*XY-XX*Y)*YZ+(XX*YY-XY*XY)*Z;

    mx = mxDET/DET;
    my = myDET/DET;
    b  = bDET/DET;

    // slopes are mx & my as is.
    // model height is given:
    float mod=mx*(p->x-x0)+my*(p->y-y0)+b+z0;
    if (1) {
      float phi=atan2(p->y,p->x);
      float mr=mx*cos(phi)+my*sin(phi);
      float mp=-mx*sin(phi)+my*cos(phi);
      //      fprintf(stderr,"WHOA %f %f %g %g %g %g %g %g %g %g %g %g\n",
      //	      p->x,p->y,mod,mx,my,mr,mp,DET,I,mxDET,myDET,bDET);
    }
    *zval = mod;    *dzdx = mx;    *dzdy = my;
  }
  return;
  // output the nearest point's parameters with zeros for slopes
 output_nearest:
  *zval=*(tl_list[0]->val);  *dzdx=0;  *dzdy=0;
  if (0) {
    fprintf(stderr,"%f %f %f %f %f %f %f\n",p->x,p->y,*(tl_list[0]->val),
	    0.0,0.0,0.0,0.0);
  }
  return;
}

int init_tlc(treelink_control *tlc) {
  float *figure;
  point *p;
  FILE *fp;
  int  arraysize=512;
  int  size=512;
  char str[size];
  float r_max;

  if ((fp=fopen(tlc->treelink_file,"r"))==NULL) {
    fprintf(stderr,"can't open input file %s\n",tlc->treelink_file);
    exit(1);
  }
  // start with initial aray size
  p=(point*)malloc(arraysize*sizeof(point));
  figure=(float*)malloc(arraysize*sizeof(float));
  int i=0;
  // read the first 2 informative lines and ignore them
  fgets(str,size,fp);  fgets(str,size,fp);
  // the meat of the tnt file:
  while (fgets(str,size,fp)) {
    if (i>=arraysize) {
      arraysize *= 2;
      if (((p=(point*)realloc(p,arraysize*sizeof(point)))==NULL) ||
	  ((figure=(float*)realloc(figure,arraysize*sizeof(float)))==NULL)) {
	fprintf(stderr,"can't reallocate. exiting..\n");
	exit(1);
      }
    }
    sscanf(str,"%f %f %f",&p[i].x,&p[i].y,&figure[i]);
    i++;
  }
  fclose(fp);
  fprintf(stderr,"read a total of %d lines\n",i);
  // datasize is datasize
  int datasize=i;
  if (((p=(point*)realloc(p,datasize*sizeof(point)))==NULL) ||
      ((figure=(float*)realloc(figure,datasize*sizeof(float)))==NULL)) {
    fprintf(stderr,"can't reallocate. exiting..\n");
    exit(1);
  }

  // data is read in. proceed to set up the tree structure. do this by 
  // setting up intermediate grids, all of which should be set up to 
  // "own" several nodes on a finer grid. The whole idea is to have 
  // overlapping ownership of the nodes.

  // first find the largest radius
  {
    int *radial_sort=(int*)malloc(datasize*sizeof(int));
    for (i=0;i<datasize;i++) radial_sort[i]=i;
    qsort_r(radial_sort,datasize,sizeof(int),compare_radii,p);
    r_max=sqrt(pow(p[radial_sort[datasize-1]].x,2)+
	       pow(p[radial_sort[datasize-1]].y,2));
  }

  // can generate up the search tree not using a fixed grid of hexagonal points
  // but a progressively finer hex grid that is centered on each of the points
  // starting from the top (7 points over aperture) thru 
  // the 2nd level (up to 7 points per subaperture)  and so on. ratios between
  // 'reach' of adjacent levels should be just over 1/sqrt(3).
  // finally each node in the figure data should probably be assigned to 
  // its 2 nearest neighbor nodes, particularly if the distance ratio is close
  // to 1.

  {
    point pt={0.0,0.0};
    build_treelink_r(&tlc->tl,&pt,r_max,tlc->lvs);
  }


  // now attach each sample to the 4 nearestnodes at the lvs-1 level
  int ntl=4;
  treelink *tl_list[ntl];
  int k;

  for (i=0;i<datasize;i++) {
    tree_nearest_links(tlc->tl,&p[i],tlc->lvs-1,tl_list,ntl);
    k=0;
    while (k<ntl) {
      if (tl_list[k]==NULL) continue;
      treelink_append(tl_list[k],&p[i]);
      // store the z (data) value. this is easier than writing a separate
      // routine with similar functionality as treelink_append()
      tl_list[k]->next_link[tl_list[k]->nelem-1]->val = &figure[i];
      k++;
    }
  }
  // and prune any branches that don't lead to data at level lvs+1
  k=tlc->lvs;
  k--;
  //  treelink_prune(tlc->tl,k,3); // need more than 3 nodes attached to penultimate
  treelink_prune(tlc->tl,k,2); // need more than 2 nodes attached to penultimate
  while (k--) { // remove nodes at level k and consolidate/remap data
    treelink_prune(tlc->tl,k,0);
  }
  // free p and figure
  //  free(p);
  free(figure);
  // set tlc as initialized.
  tlc->inited=1;

  return(0);
  // now to test.
  FILE *OOP;
  if ((OOP=fopen("tmp.tnt","w"))==NULL) {
    fprintf(stderr,"can't open file for writing.\n exiting..\n");
    exit(1);
  }

  for (i=0;i<datasize;i++) {
    int nn=tree_nearest_links(tlc->tl,&p[i],tlc->lvs,tl_list,ntl);
    if (nn>2) { 
      /* should normally be the case, since tree was pruned >2 above */
      /* can do vector calculations maybe. */
      /* with model z(x,y) = mx*x + my*y + b */
      
      /* /x2 xy x\ /mx\   / xz \ */
      /* |xy y2 y| |my| = | yz | */
      /* \x  y  I/ \b /   \ z  / */
      /* inverse of the sampling matrix is the transpose of cofactor matrix */
      /* divided through by the determinant of the sampling matrix. */
      float XX=0,XY=0,YY=0,X=0,Y=0,I=0,XZ=0,YZ=0,Z=0;
      float x,y,z,x0,y0,z0;
      x0=tl_list[0]->pt->x;
      y0=tl_list[0]->pt->y;
      z0=*(tl_list[0]->val);
      for (k=0;k<nn;k++) {
	x = tl_list[k]->pt->x - x0;
	y = tl_list[k]->pt->y - y0;
	z = *(tl_list[k]->val)- z0;
	XX += x*x;	XY += x*y;	YY += y*y;
	X  += x;	Y  += y;	I  += 1;
	XZ += x*z;	YZ += y*z;	Z  += z;
      }
      float DET=XX*(YY*I-Y*Y)+XY*(Y*X-I*XY)+X*(XY*Y-X*YY);
      if (DET==0) {
	//	fprintf(stderr,"youch! DET is zero..\n");
	goto output_nearest;
      }
      float mx,my,b;
      mx=(YY*I-Y*Y)*XZ/DET+(X*Y-I*XY)*YZ/DET+(XY*Y-X*YY)*Z/DET;
      my=(Y*X-I*XY)*XZ/DET+(I*XX-X*X)*YZ/DET+(X*XY-XX*Y)*Z/DET;
      b =(XY*Y-YY*X)*XZ/DET+(X*XY-XX*Y)*YZ/DET+(XX*YY-XY*XY)*Z/DET;
      // slopes are mx & my as is.
      // model height is given:
      float mod=mx*(p[i].x-x0)+my*(p[i].y-y0)+b+z0;
      float phi=atan2(p[i].y,p[i].x);
      float mr=mx*cos(phi)+my*sin(phi);
      float mp=-mx*sin(phi)+my*cos(phi);
      fprintf(OOP,"%f %f %g %g %g %g %g\n",p[i].x,p[i].y,mod,mx,my,mr,mp);
    } 
    continue;
    // output the nearest point's parameters with zeros for slopes
  output_nearest:
    fprintf(OOP,"%f %f %f %f %f %f %f\n",p[i].x,p[i].y,*(tl_list[i]->val),
	    0.0,0.0,0.0,0.0);

  }
  fclose(OOP);
  // and print nodes all the way out to data
  //  int m;
  //  for (m=0;m<lvs+1;m++) tree_printlevel(tl_head,m);
  return(0);
}

void build_treelink_r(treelink **tl,point *p,float radius,int lvl) {
  // recursive routine to populate 2D tree
  int i;
  point pt2;
  treelink_init(tl,p);
  if (lvl==0) return;
  for (i=0;i<6;i++) {
    pt2.x=p->x+radius*cos(i*M_PI/3.);
    pt2.y=p->y+radius*sin(i*M_PI/3.);
    treelink_append(*tl,&pt2);
    build_treelink_r(&((*tl)->next_link[i]),&pt2,radius/sqrt(3),lvl-1);
  }
}

point *atree_interp(treelink *tl,point *p) {
  // returns nearest node on the end of the treelink to point p.
  // do this recursively until no more nodes are available.
  int i,min_i;
  float d2min,d2;
  if (tl->nelem>0) {
    d2min=pow(p->x-tl->points[0].x,2)+pow(p->y-tl->points[0].y,2);
    min_i=0;
    for (i=1;i<tl->nelem;i++) {
      d2=pow(p->x-tl->points[i].x,2)+pow(p->y-tl->points[i].y,2);
      if (d2<d2min) {
	d2min=d2;
	min_i=i;
      }
    }
    // check for daughter nodes at min_i. return this point or continue down
    if (tl->next_link[min_i]->nelem==0) {
      return(&(tl->points[min_i]));
    } else {
      return(atree_interp(tl->next_link[min_i],p));
    }
  } else {
    fprintf(stderr,"gave me a node with no petals.\n");
    return(NULL);
  }
}

void tree_printlevel(treelink *tl,int lvl) {
  int j;
  static int mylvl=0;
  if (lvl==0) {
    for (j=0;j<tl->nelem;j++) {
      fprintf(stderr,"%d %f %f\n",mylvl,tl->points[j].x,tl->points[j].y);
    }
    return;
  } else {
    for (j=0;j<tl->nelem;j++) {
      mylvl++;
      tree_printlevel(tl->next_link[j],lvl-1);
      mylvl--;
    }
  }
}

treelink *atree_nearest_link(treelink *tl,point *p,int lvl) {
  // returns nearest link node on the end of the treelink to point p.
  // do this only down to lvl from current.
  int i,min_i;
  float d2min,d2;
  if (tl->nelem>0) {
    d2min=pow(p->x-tl->points[0].x,2)+pow(p->y-tl->points[0].y,2);
    min_i=0;
    for (i=1;i<tl->nelem;i++) {
      d2=pow(p->x-tl->points[i].x,2)+pow(p->y-tl->points[i].y,2);
      if (d2<d2min) {
	d2min=d2;
	min_i=i;
      }
    }
    // check for daughter nodes at min_i. return this point or continue down
    if ((lvl==0) || (tl->next_link[min_i]->nelem==0)) {
      return(tl->next_link[min_i]);
    } else {
      return(atree_nearest_link(tl->next_link[min_i],p,lvl-1));
    }
  } else {
    fprintf(stderr,"gave me a node with no petals.\n");
    return(NULL);
  }
}

int tree_nearest_links(treelink *tl,point *p,int lvl,treelink **tl_list,int ntl) {
  // populate tl_list[0..ntl-1] with the treelink*s nearest to point *p
  // at lvl. done recursively.

  int i,min_i;
  float d2min,d2;
  if (lvl>0) {
    if (tl->nelem>0) {
      d2min=pow(p->x-tl->points[0].x,2)+pow(p->y-tl->points[0].y,2);
      min_i=0;
      for (i=1;i<tl->nelem;i++) {
	d2=pow(p->x-tl->points[i].x,2)+pow(p->y-tl->points[i].y,2);
	if (d2<d2min) {
	  d2min=d2;
	  min_i=i;
	}
      }
      return(tree_nearest_links(tl->next_link[min_i],p,lvl-1,tl_list,ntl));
    } else {
      fprintf(stderr,"gave me a node with no petals.\n");
      return(0);
    }
  } else {
    // if ntl > tl->nelem fill extras with NULL
    int ixlist[tl->nelem];
    float dlist[tl->nelem];

    if (ixlist==NULL || dlist==NULL) {
      fprintf(stderr,"can't allocate!! exiting..\n");
      exit(1);
    }
    int k;
    for (k=0;k<tl->nelem;k++) {
      ixlist[k]=k;
      dlist[k] = pow(p->x-tl->points[k].x,2)+pow(p->y-tl->points[k].y,2);
    }
    qsort_r(ixlist,tl->nelem,sizeof(int),compare_float_array,dlist);
    for (k=0;k<ntl;k++) {
      if (k<tl->nelem) {
    	tl_list[k]=tl->next_link[ixlist[k]];
      } else {
    	tl_list[k]=NULL;
      }
    }
    return(((ntl<tl->nelem)?ntl:tl->nelem));
  }
}

void treelink_init (treelink **tl,point *p) {
  // allocate
  *tl = (treelink*)malloc(sizeof(treelink));
  // copy over point value
  if (p==NULL) {
    (*tl)->pt=NULL;
  } else {
    (*tl)->pt=(point*)malloc(sizeof(point));
    memcpy((*tl)->pt,p,sizeof(point));
  }

  (*tl)->val=NULL;
  (*tl)->maxelem=64;
  (*tl)->nelem=0;
  (*tl)->points=(point*)malloc(((*tl)->maxelem)*sizeof(point));
  (*tl)->next_link=(treelink**)malloc(((*tl)->maxelem)*sizeof(treelink*));
  int i;
  for (i=0;i<(*tl)->maxelem;i++) {
    (*tl)->next_link[i]=(treelink*)NULL;
  }
}

void treelink_append (treelink *tl,point *p) {
  // if necessary, reallocate
  if (tl->nelem >= tl->maxelem) {
    tl->maxelem *= 2;
    tl->points=(point*)realloc(tl->points,tl->maxelem*sizeof(point));
    tl->next_link=(treelink**)realloc(tl->next_link,
				      tl->maxelem*sizeof(treelink*));
  }
  // and continue append
  memcpy(&(tl->points[tl->nelem]),p,sizeof(point));
  treelink_init(&(tl->next_link[tl->nelem]),p);
  tl->nelem++;
}

void treelink_prune (treelink *tl,int lvl,int threshold) {
  // work down to level lvl
  int i;
  if (lvl==0) {
    // starting here, inspect contents of daughter nodes and trim as necessary
    i=tl->nelem;
    while (i--) {
      if (tl->next_link[i]->nelem<=threshold) {
	// remove tl->next_link[i] and tl->points[i]
	// but shift the array contents above it to start
	treelink_destroy(tl->next_link[i]);
	if (i<tl->nelem-1) {
	  memmove(&tl->next_link[i],
		  &tl->next_link[i+1],(tl->nelem-i-1)*sizeof(treelink*));
	  memmove(&tl->points[i],
		  &tl->points[i+1],(tl->nelem-i-1)*sizeof(point*));
	}
	tl->nelem--;
      }
    }
  } else {
    for (i=0;i<tl->nelem;i++)
      treelink_prune(tl->next_link[i],lvl-1,threshold);
  }
}

void treelink_destroy(treelink *tl) {
  // reverse allocations in treelink_init for tl
  // any non-deallocated daughter nodes may be left dangling
  int i;
  for (i=0;i<tl->nelem;i++) treelink_destroy(tl->next_link[i]);
  free(tl->next_link);
  free(tl->points);
  free(tl->pt);
  free(tl);
}

int compare_float_array(const void* a,const void* b,void* arg) {
  float *fa=(float*)arg;
  int   ia=*(int*)a;
  int   ib=*(int*)b;
  if (fa[ia]<fa[ib]) return(-1);
  if (fa[ia]>fa[ib]) return(+1);
  return(0);
}

int compare_radii(const void* a,const void* b,void* arg) {
  float r2a,r2b;
  point *p=(point*)arg;
  r2a = pow(p[*(int*)a].x,2)+pow(p[*(int*)a].y,2);
  r2b = pow(p[*(int*)b].x,2)+pow(p[*(int*)b].y,2);
  if (r2a < r2b) return(-1);
  if (r2a > r2b) return(1);
  return(0);
}
