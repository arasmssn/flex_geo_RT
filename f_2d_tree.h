typedef struct {
  float x,y;
} point;

typedef struct treelink {
  point           *pt;
  float          *val;
  int             maxelem;
  int             nelem;
  point           *points;
  struct treelink **next_link;
} treelink;

typedef struct {
  char *treelink_file;
  int  inited;
  int  lvs;
  treelink *tl;
} treelink_control;
 
void tl_interp(treelink_control *tlc,point *p,float *z, float *dzdx, float *dzdy);
int  init_tlc(treelink_control *tlc);
void treelink_destroy(treelink *tl);
void treelink_prune (treelink *tl,int lvl,int threshold);
void tree_printlevel(treelink *tl,int lvl);
treelink *atree_nearest_link(treelink *tl,point *p,int lvl);
int  tree_nearest_links(treelink *tl,point *p,int lvl,treelink **tl_list,int ntl);
point    *atree_interp(treelink *tl,point *p);
void build_treelink_r(treelink **tl,point *p,float radius,int lvl);
void treelink_init (treelink **tl,point* p);
void treelink_append (treelink *tl,point *p);
int  compare_float_array(const void* a,const void* b,void* arg);
int  compare_radii(const void* a,const void* b,void* arg);
