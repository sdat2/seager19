/* $Id: tios.h,v 1.2 1997/03/21 20:37:32 senya Exp senya $
   $Revision: 1.2 $
   $Date: 1997/03/21 20:37:32 $
   $Author: senya $
*/
/********* header file for TIOS i/o-system ***************/
#define Magic "TiOs!9(1"
#define Version 11
#define TIOS_RCSID = "$Id: tios.h,v 1.2 1997/03/21 20:37:32 senya Exp senya $"

#define ABS(x) ((x)>=0?(x):-(x))
#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define APPROX_EQ(x,y,e) ((float)ABS((x)-(y)) < 0.5*e)
#define ENSO_DAY (12./365.)

#define REAL float
#define UINT unsigned int

#define FRST_WORD "first"
#define LAST_WORD "last"
#define AUTO_FLAG (-987654)

#define MAX_SLICES 50

#define TP_STREAM  0
#define TP_CLUSTER 1

#define OLD_TIOS_HEADER 11
#define NEW_TIOS_HEADER 16

#define FSKIP(f,s) fseek(f,(size_t)(s),SEEK_CUR)

typedef struct ADDR
{
  UINT val;
  REAL time;
  struct ADDR *next;
} ADDR;

typedef struct RG_ST
{
  int num, fmt, absolute;
  REAL sta, end, step;
}RG_ST;

typedef struct CLUSTER
{
  int   type, n, nid, ntag;    /* type, size of a cluster, IDs, Tags */
  REAL *px, *py;               /* grids */
  UINT *pid;                   /* IDs (optional) */
  char *ptag;                  /* Tags (optional) */
  UINT *indx;                  /* indexes of data */
}CLUSTER;

typedef struct RANGE
{
  struct RG_ST x, y, z, t;
  int    range_nu, type;
  struct CLUSTER cluster;
  struct RANGE *next;
}RANGE;

typedef struct MAP
{
  int map_nu, type;
  int NMAP, NXY;
  int *cmp;
  struct MAP *next;
}MAP;

typedef struct GRID
{
  int grid_nu, type;
  int NX, NY, NZ;
  REAL *xx, *yy, *zz;
  struct GRID *next;
}GRID;

typedef struct GR_MAP
{
  int gmap_nu, type;
  struct GRID  *grid;
  struct RANGE *range;
  int MX, MY, MZ;
  int *mapx, *mapy, *mapz;
  struct GR_MAP *next;
}GR_MAP;

#define LABEL_STREAM_LEN 48
#define MAX_ADJS 32

typedef struct STREAM
{
  int stream_nu, type;
  char label[LABEL_STREAM_LEN];
  struct {
    char id;                           
    int  num, *adr;
    REAL *val;
  } slice;
  struct {
    int count, init, fmt, absolute;    /* number of time points; time flags */
    REAL sta, end, step, next, *data;  /* time grid parameters */
  } time;
  int size;
  struct GR_MAP *gmap;
  struct {
    int count, num[MAX_ADJS];     /* number of vars in a stream; IDs of vars */
    struct VAR *adr[MAX_ADJS];    /* pointers to VARS */
  } vars;
  struct CLUSTER cluster;
  UINT *laddr;
  struct ADDR *addr, *curr;
  struct STREAM *next;
} STREAM;

#define STR_VAR_LEN 12
typedef struct VAR
{
  int  var_nu, type, blength, used;   
  char label[STR_VAR_LEN];
  REAL *var;                           /* data pointer */
  REAL flag;
  void (*func)();                      /* virtual function */
  struct MAP  *map;
  struct GRID *grid;
  struct {                             /* for a base-var: list of streams */ 
    int count, num[MAX_ADJS];
    struct STREAM *adr[MAX_ADJS];
  } strs;
  struct VAR *next;
} VAR;

#define TITLE_TIOS_LEN 200

typedef struct TIOS 
{
  int version;
  int map_count;
  int grid_count;
  int gmap_count;
  int rang_count;
  int str_count;
  int clust_count;
  int var_count;
  int time_fmt;
  REAL pole_alpha,pole_beta,pole_gamma;
  struct {
    int everystep;
    int debug;
    int updated;
    int restart;
  } cntrl;
  REAL time_begin, time_end, time_resolution;
  struct MAP      *maps;
  struct GRID     *grids;
  struct GR_MAP   *gr_maps;
  struct RANGE    *ranges;
  struct VAR      *vars;
  struct STREAM   *streams;
  char tios_name[TITLE_TIOS_LEN];
  long addr_start;
  long addr_grids;
  long addr_vars;
  long addr_strs;
  long addr_end;
}TIOS;

#ifdef TIOS_APPL

VAR    *var_curr;
GRID   *grid_curr;
GR_MAP *gmap_curr;
STREAM *str_curr;
ADDR   *addr_curr;
RANGE  *rang_curr;
MAP    *map_curr;

FILE   *indx_file, *data_file, *tios_file, *ingr_file;

TIOS   tios;

#else

static VAR    *var_curr;
static GRID   *grid_curr;
static GR_MAP *gmap_curr;
static STREAM *str_curr;
static ADDR   *addr_curr;
static RANGE  *rang_curr;
static MAP    *map_curr;

static FILE   *indx_file, *data_file, *tios_file, *ingr_file;

static TIOS   tios;

#endif
/********************************************/









