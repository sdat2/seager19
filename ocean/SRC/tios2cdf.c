/****************************************************************************
*  This is a TIOS -> NetCDF conversion program. 
****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <values.h>
#include <ctype.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdint.h>

#include "tios.h"
#include "cuf.h"

#define F_INDX "tios.indx"
#define F_DATA "tios.data"
#define F_OUT  "tios.nc"
#define F_NAME  "tios"

static int   if_table, if_debug, if_old, if_netcdf, if_netcdf4grads;
static int   if_list, Pick;
static char *f_indx, *f_data, *f_out, *f_name;
static char *DumpStream, *DumpVar, *DumpFullVar;
int  nbegin, nend, njump;
float  fmiss = AUTO_FLAG;

#define IF_TABLE if(if_table)
#define IF_LIST if(if_list)

static int flen(str)
char *str;
{
  char *pstr = str;
  while 
    (
     (int)*pstr           && 
     !iscntrl((int)*pstr) && 
     !isspace((int)*pstr)
    ) pstr++;

  return (int)(pstr-str);
}

static int *real_to_int(n, real)
int n;
REAL *real;
{
  REAL *pr;
  int *pi, *pp;
  
  pi = pp = (int *) malloc (n*sizeof(int));
  for (pr = real; pr < real + n; pr++) *pp++ = (int)(*pr);
  return pi;
}

static REAL *fread_alloc (n, file)
int n;
FILE *file;
{
  REAL *data;

  if (n) 
    {
      data = (REAL *) malloc(n*sizeof(REAL));
      loc_rd (sizeof(REAL), n, data, file);
      return data;
    }
  else
    return (REAL *)NULL;
}

static int *rd_as_real (n, file)
int n;
FILE *file;
{
  REAL *rdata;
  int  i, *idata;

  if (n) 
    {
      idata = (int *) malloc(n*sizeof(int));
      rdata = (REAL *) malloc(n*sizeof(REAL));
      loc_rd (sizeof(REAL), n, rdata, file);
      for (i = 0; i < n; i++) idata[i] = rdata[i];
      free (rdata);
      return idata;
    }
  else
    return (int *)NULL;
}

static VAR *var_by_number (n)
int n;
{
  VAR *pvar;

  pvar = tios.vars;
  while (pvar->next && (pvar->var_nu != n) ) pvar = pvar->next;
  return pvar;
}


static void lread_tios (file)
FILE *file;
{
  static REAL tmp[NEW_TIOS_HEADER];
  static char magic[8];

  loc_rd (1, 8, magic, file);
  if ( ! strncmp(magic, Magic, 8) ) { /* new version */
    loc_rd (sizeof(REAL), 1, (void *)tmp, file);
    tios.version = (int)tmp[0];
  }
  else {                              /* an old 0-version */
    fseek (file, 0, SEEK_SET);
    tios.version = 0;
  }

  if ( tios.version == 0) 
    loc_rd (sizeof(REAL), OLD_TIOS_HEADER,  (void *)tmp, file);
  else {
    loc_rd (sizeof(REAL), NEW_TIOS_HEADER,  (void *)tmp, file);
    tios.clust_count  = (int)tmp[11];
    tios.pole_alpha  = tmp[12];
    tios.pole_beta  =  tmp[13];
    tios.pole_gamma  = tmp[14];
  }

  tios.grid_count = (int)tmp[0];
  tios.gmap_count = (int)tmp[1];
  tios.var_count  = (int)tmp[2];
  tios.str_count  = (int)tmp[3];
  tios.time_begin =      tmp[4];
  tios.time_end   =      tmp[5];
  tios.addr_start = (int)tmp[6];
  tios.addr_grids = (int)tmp[7];
  tios.addr_vars  = (int)tmp[8];
  tios.addr_strs  = (int)tmp[9];
  tios.addr_end   = (int)tmp[10];

  fread (tios.tios_name, (size_t)1, (size_t)TITLE_TIOS_LEN, file);

  IF_TABLE 
    {
      printf ("TIOS - Transfer IO System; Version %d.%d\n\n", 
	      tios.version/10, tios.version%10);
      printf ("Label     : %s\n", tios.tios_name);
      printf ("Grids     : %d, (mapped %d).\n",tios.grid_count,tios.gmap_count);
      printf ("Variables : %d\n", tios.var_count);
      printf ("Streams   : %d\n", tios.str_count);
      if (tios.clust_count) printf ("Clusters   : %d\n", tios.clust_count);
      printf ("Start time : %f; End time : %f\n",tios.time_begin,tios.time_end);
      printf ("Pole_alpha : %f\n", tios.pole_alpha);
      printf ("Pole_beta  : %f\n", tios.pole_beta);
      printf ("Pole_gamma : %f\n", tios.pole_gamma);
    }
  fseek (file, (size_t)tios.addr_grids, SEEK_SET);
  tios.grids   = grid_curr = (GRID *)calloc (1, sizeof(GRID)); 
  tios.gr_maps = gmap_curr = (GR_MAP *)calloc (1, sizeof(GR_MAP)); 
  tios.vars    = var_curr  = (VAR *)calloc (1, sizeof(VAR)); 
  tios.streams = str_curr  = (STREAM *)calloc (1, sizeof(STREAM)); 
}

static void lread_grid (ng)
int ng;
{
  static REAL tmp[3];
  REAL *xx, *yy, *zz, *real;

  loc_rd (sizeof(REAL), 3, tmp, indx_file);

  xx = fread_alloc ((int)tmp[0], indx_file);
  yy = fread_alloc ((int)tmp[1], indx_file);
  zz = fread_alloc ((int)tmp[2], indx_file);

  IF_TABLE
    {
      printf ("Grid #%d: -------------------------------\n", ng);
      printf (" NX = %3d", (int)tmp[0]);
      if ((int)tmp[0]) 
	printf ("\txx : %.4g %.4g %.4g ... %.4g\n", 
		xx[0],xx[1],xx[2],xx[(int)tmp[0]-1]);
      printf (" NY = %3d", (int)tmp[1]);
      if ((int)tmp[1]) 
	printf ("\tyy : %.4g %.4g %.4g ... %.4g\n", 
		yy[0],yy[1],yy[2],yy[(int)tmp[1]-1]);
      printf (" NZ = %3d", (int)tmp[2]);
      if ((int)tmp[2]) 
	printf ("\tzz : %.4g %.4g %.4g ... %.4g\n", 
		zz[0],zz[1],zz[2],zz[(int)tmp[2]-1]);
    }

  grid_curr->grid_nu = ng;
  grid_curr->NX = (int)tmp[0];
  grid_curr->NY = (int)tmp[1];
  grid_curr->NZ = (int)tmp[2];

  grid_curr->xx = xx;
  grid_curr->yy = yy;
  grid_curr->zz = zz;

  grid_curr = (grid_curr->next = (GRID *)calloc (1, sizeof (GRID)));
}

static void lread_gmap (nmap)
int nmap;
{
  GRID *pgr;
  static REAL tmp[4];

  loc_rd (sizeof(REAL), 4, tmp, indx_file);
  
  gmap_curr->gmap_nu = nmap;
  gmap_curr->MX = (int)tmp[0];  gmap_curr->MY = (int)tmp[1];
  gmap_curr->MZ = (int)tmp[2];  

  pgr = tios.grids;
  while (!pgr->next && pgr->grid_nu != (int)tmp[3]) pgr = pgr->next;
  if (!pgr->next) {
    printf ("TIOS: inconsistent information for grid map# %d.Stop.\n", nmap);
    exit (-1);
  }
  gmap_curr->grid = pgr;

  if (gmap_curr->MX) gmap_curr->mapx = rd_as_real(gmap_curr->MX, indx_file);
  if (gmap_curr->MY) gmap_curr->mapy = rd_as_real(gmap_curr->MY, indx_file);
  if (gmap_curr->MZ) gmap_curr->mapz = rd_as_real(gmap_curr->MZ, indx_file);

  IF_TABLE
    printf ("Grid Map #%d: MX = %3d, MY = %3d, MZ = %3d\n", 
	    nmap, gmap_curr->MX, gmap_curr->MY, gmap_curr->MZ);

  gmap_curr = (gmap_curr->next = (GR_MAP *)calloc (1, sizeof (GR_MAP)));
}

static void lread_var (nvar)
int nvar;
{
  int i;
  static REAL tmp[3], *pre;
  char label[STR_VAR_LEN];

  var_curr->var_nu = nvar;
  loc_rd (sizeof(REAL), 3, tmp, indx_file);
  var_curr->strs.count = (int)tmp[0];

  if ( var_curr->strs.count ) 
    {
      pre = (REAL *) malloc (var_curr->strs.count * sizeof(REAL));
      loc_rd (sizeof(REAL), var_curr->strs.count, pre, indx_file);
      for (i = 0; i < var_curr->strs.count; i++) var_curr->strs.num[i] = pre[i];
      free (pre);
    }
    
  fread (label, (size_t)1, (size_t)STR_VAR_LEN, indx_file);
  IF_TABLE {
    printf ("Variable #%d <%s>\tis the Base for %d streams:", 
	    nvar, label, (int)tmp[0]);
    for (i = 0; i < var_curr->strs.count; i++) 
      printf (" %d", var_curr->strs.num[i]);
    putchar ('\n');
  }
  strcpy (var_curr->label, label);

  var_curr->grid = tios.grids;
  while ((int)tmp[1] != var_curr->grid->grid_nu) 
    var_curr->grid = var_curr->grid->next;

  var_curr->flag = tmp[2];
  var_curr = (var_curr->next = (VAR *)calloc (1, sizeof(VAR)));
}

static void lread_stream (nstream)
int nstream;
{
  int i, ngmap, *pint;
  GR_MAP *pgm;
  VAR    *pvar;
  static REAL tmp[9];

  if (tios.version) {
    loc_rd (sizeof(REAL), 1, tmp, indx_file);
    str_curr->type = (int)tmp[0];
  }
  else
    str_curr->type = 0;

  str_curr->stream_nu = nstream;
  loc_rd (sizeof(REAL), 9, tmp, indx_file);

  str_curr->time.count = (int)tmp[0];
  str_curr->time.sta   =      tmp[1];
  str_curr->time.end   =      tmp[2];
  str_curr->time.step  =      tmp[3];

  str_curr->slice.id   =(char)tmp[4];
  str_curr->slice.num  = (int)tmp[5];
  ngmap                = (int)tmp[6]; 
  str_curr->size       = (int)tmp[7]; 
  str_curr->vars.count = (int)tmp[8];

  str_curr->slice.val = 
    (REAL *)malloc (str_curr->slice.num * sizeof(REAL));
  loc_rd (sizeof(REAL), str_curr->slice.num, str_curr->slice.val, indx_file);
  
  pint = rd_as_real (str_curr->vars.count, indx_file);
  for (i = 0; i < str_curr->vars.count; i++) 
    str_curr->vars.adr[i] = var_by_number(str_curr->vars.num[i] = pint[i]);  
  free (pint);

  fread (str_curr->label, (size_t)1, (size_t)LABEL_STREAM_LEN, indx_file);

  IF_TABLE {
    printf ("Stream%s #%d: label:<%s>; records = %d\n %c = [ ",
      ( str_curr->type == TP_CLUSTER ? "/Cluster" : ""),
      nstream, str_curr->label, str_curr->time.count, str_curr->slice.id);
    for (i = 0; i < str_curr->slice.num; i++)
      printf ("%.2f ", str_curr->slice.val[i]);
    printf ("]\n vars = %d [", str_curr->vars.count); 
    for (i = 0; i < str_curr->vars.count; i++) 
      printf (" %s ", str_curr->vars.adr[i]->label);
    if ( str_curr->type == TP_CLUSTER ) 
      printf ("] npts=%d\n", str_curr->cluster.n);
    else printf ("]\n");
  }
  IF_LIST {
    printf ("tios2cdf -s %s -o STREAM_%s.nc -f %s\n", 
	    str_curr->label, str_curr->label, f_name);
  }

  pgm = tios.gr_maps;
  
  while (pgm->next && ngmap != pgm->gmap_nu) pgm = pgm->next;
  if (!ngmap || !pgm->next) 
    str_curr->gmap = NULL;
  else
    str_curr->gmap = pgm;
  
  str_curr->laddr     = (UINT *) malloc(sizeof(UINT) * str_curr->time.count);
  str_curr->time.data = (REAL *) malloc(sizeof(REAL) * str_curr->time.count);

  if (if_old) {
    for (i = 0; i < str_curr->time.count; i++) {
      loc_rd (sizeof(REAL), 2, tmp, indx_file);
      str_curr->laddr[i]     = (UINT)tmp[0];
      str_curr->time.data[i] = tmp[1];
    }
  }
  else {
    loc_rd (sizeof(UINT), str_curr->time.count, str_curr->laddr, indx_file);
    loc_rd (sizeof(REAL), str_curr->time.count, str_curr->time.data, indx_file);
  }

  str_curr->next = (STREAM *)calloc (1, sizeof (STREAM));
}

static void do_table ()
{
  int ngr, nvar, nstream, nadr, ladr;
  
  lread_tios (indx_file);

  for (ngr = 1; ngr <= tios.grid_count; ngr++) lread_grid (ngr);
  IF_TABLE putchar ('\n');

  for (ngr = 1; ngr <= tios.gmap_count; ngr++) lread_gmap (ngr);
  IF_TABLE putchar ('\n');

  for (nvar = 1; nvar <= tios.var_count; nvar++) lread_var (nvar);
  IF_TABLE putchar ('\n');

  for (nstream = 1; nstream <= tios.str_count; nstream++) {
    lread_stream(nstream);
    str_curr = str_curr->next;
  }
  IF_TABLE 
  printf ("----------------------------End of Contents for file <%s>\n",f_data);
}

static REAL *remap_grid (n, map, x)
int n, *map;
REAL *x;
{
  int i;
  REAL *tmp;
  
  tmp = (REAL *) malloc (n * sizeof (REAL));
  for (i = 0; i <  n; i++) tmp[i] = x[map[i]];

  return tmp;
}

#define STR77(s) s,(long)strlen(s)

int do_cdf ()
{
  int idf, i = 2, ONE = 1;
  STREAM *pstr;

    odb_open_(&idf, f_out, &i, strlen(f_out)); 

    odb_setcattr_(&idf, "global", "model", f_data, 
		  6, 5, strlen(f_data));
    odb_setrattr_(&idf, "global", "missing_value", &fmiss, 6, 13);

    odb_dfgr_(&idf, "GRID", &ONE, 4);
    odb_setrattr_(&idf, "GRID", "pole_alpha", &tios.pole_alpha, 4, 10);
    odb_setrattr_(&idf, "GRID", "pole_beta",  &tios.pole_beta, 4, 9);
    odb_setrattr_(&idf, "GRID", "pole_gamma", &tios.pole_gamma, 4, 9);

  if (if_debug) printf ("Output File: %s\n", f_out);
  pstr = tios.streams;
  while (pstr->next) {
    
    if ( ! pstr->time.count ) {pstr = pstr->next; continue;}
 
    if ( DumpStream ) {
      if (! strcmp (DumpStream, pstr->label)) {
	do_stream (idf, pstr);
	break;
      }
    }
    else do_stream (idf, pstr);

    pstr = pstr->next;
  }
  
    odb_close_(&idf);
}

void def_stream (idf, pstr, pcnt, s_tag, v_tag, t_tag)
int     idf;
STREAM *pstr;
char   *pcnt, *s_tag, *v_tag, *t_tag;
{
  VAR    *pvar;
  GRID   *pgr;
  GR_MAP *pgm;
  char    xyzt[5], *p1;
  REAL   *xx, *yy, *zz;
  int     i, NX, NY, NZ, def_grids = 0;

  pgm = pstr->gmap;
  pgr = pgm ? pgm->grid : pstr->vars.adr[0]->grid;

  /* determine grids parameters: */
  switch (pstr->slice.id) {
  case 'X': case 'x':
/*    if (DumpVar)
       strcpy (xyzt, "YZXT");
    else*/
       strcpy (xyzt, "YLXT");
    NX = pstr->slice.num; xx = pstr->slice.val;
      if (pgm && pgm->MY) 
	{ NY = pgm->MY; yy = remap_grid (pgm->MY, pgm->mapy, pgr->yy); }
      else 
	{ NY = pgr->NY; yy = pgr->yy; }
    if (pgm && pgm->MZ) 
      { NZ = pgm->MZ; zz = remap_grid (pgm->MZ, pgm->mapz, pgr->zz); }
    else 
      { NZ = pgr->NZ; zz = pgr->zz; }
    break;
  case 'Y': case 'y':
/*    if (DumpVar)
       strcpy (xyzt, "XZYT");
    else*/
       strcpy (xyzt, "XLYT");
    NY = pstr->slice.num; yy = pstr->slice.val;
    if (pgm && pgm->MX) 
      { NX = pgm->MX; xx = remap_grid (pgm->MX, pgm->mapx, pgr->xx); }
    else 
      { NX = pgr->NX; xx = pgr->xx; }
    if (pgm && pgm->MZ) 
      { NZ = pgm->MZ; zz = remap_grid (pgm->MZ, pgm->mapz, pgr->zz); }
    else 
      { NZ = pgr->NZ; zz = pgr->zz; }
    break;
  case 'Z': case 'z':
  case 'L': case 'l':
/*    if (DumpVar)
       strcpy (xyzt, "XYZT");
    else*/
       strcpy (xyzt, "XYLT");
    NZ = pstr->slice.num; zz = pstr->slice.val;
    if (pgm && pgm->MX) 
      { NX = pgm->MX; xx = remap_grid (pgm->MX, pgm->mapx, pgr->xx); }
    else 
      { NX = pgr->NX; xx = pgr->xx; }
    if (pgm && pgm->MY) 
      { NY = pgm->MY; yy = remap_grid (pgm->MY, pgm->mapy, pgr->yy); }
    else 
      { NY = pgr->NY; yy = pgr->yy; }
    break;
  }

  /* define all stream's variables */
  for (i = 0; i < pstr->vars.count; i++) { 
    pvar = pstr->vars.adr[i];

    v_tag = strcat(strcpy(v_tag, s_tag), pvar->label);
    t_tag = pvar->label;
    
    if ( DumpVar ) {
      if ( !strcmp (DumpVar, pvar->label) ) {
	if ( !def_grids ) do_grids (idf, pcnt, NX, xx, NY, yy, NZ, zz, 
				  pstr->time.count, pstr->time.data);
	def_grids = 1;
	if (if_debug) printf ("Set Variable: <%s>\n", t_tag);
	def_var(idf, t_tag, v_tag, (intptr_t)pgm, pcnt, xyzt); 
	break;
      }
    }
    else {if ( DumpFullVar ) {
      if (if_debug) printf ("compare variables: <%s> <%s>\n", DumpFullVar,v_tag);
      if ( !strcmp (DumpFullVar, v_tag) ) {
	if ( !def_grids ) do_grids (idf, pcnt, NX, xx, NY, yy, NZ, zz, 
				  pstr->time.count, pstr->time.data);
	def_grids = 1;
	if (if_debug) printf ("Set Variable: <%s> <%s>\n", t_tag,v_tag);
	def_var(idf, t_tag, v_tag, (intptr_t)pgm, pcnt, xyzt); 
	break;
      }
      }
    else { if ( DumpStream ) {
      if (!def_grids) do_grids (idf, pcnt, NX, xx, NY, yy, NZ, zz, 
				pstr->time.count, pstr->time.data);
      def_grids = 1;
      if (if_debug) printf ("Set Variable: <%s>\n", v_tag);
      def_var(idf, t_tag, t_tag, (intptr_t)pgm, pcnt, xyzt);
      }
    else {
      if (!def_grids) do_grids (idf, pcnt, NX, xx, NY, yy, NZ, zz, 
				pstr->time.count, pstr->time.data);
      def_grids = 1;
      if (if_debug) printf ("Set Variable: <%s>\n", v_tag);
      def_var(idf, v_tag, t_tag, (intptr_t)pgm, pcnt, xyzt);
    }}
  }
  }
}

void write_stream (idf, pstr, pcnt, s_tag, v_tag)
int     idf;
STREAM *pstr;
char   *pcnt, *s_tag, *v_tag;
{
  char buf[16];
  VAR    *pvar;
  REAL   *pbuf;
  int     i, j, bsize, stream = (pstr->type == TP_STREAM);

  bsize = pstr->slice.num * (stream ? pstr->size : pstr->cluster.n);
  pbuf  = (REAL *) malloc (bsize * sizeof(REAL));
  
  for (j = 0; j < pstr->time.count; j++) {
    fseek (data_file, pstr->laddr[j], SEEK_SET);
    
    for (i = 0; i < pstr->vars.count; i++) {
      pvar = pstr->vars.adr[i];      

      /* make Variable name concatenating: "STREAM.VAR" */

      if ( DumpVar ) {
        v_tag = pvar->label;
	if ( !strcmp(DumpVar, pvar->label) ) {
	  if (if_debug) printf ("Write Var: <%s> indx= %d\n", v_tag, j+1);
          rdwr_var (idf, j+1, bsize, v_tag, pvar, pbuf);
          sprintf (buf, "T\0", pcnt);
          odb_wrgr_(&idf, buf , pstr->time.data , strlen(buf));
	  break;
	}
	else
	  fseek(data_file, pstr->size*pstr->slice.num*sizeof(REAL), SEEK_CUR);
      }
      else {
      if ( DumpFullVar ) {
        v_tag = strcat(strcpy(v_tag, s_tag), pvar->label);
	if ( !strcmp(DumpFullVar, v_tag) ) {
	  if (if_debug) printf ("Write Var: <%s> indx= %d\n", v_tag, j+1);
          rdwr_var (idf, j+1, bsize, v_tag, pvar, pbuf);
	  break;
	}
	else
	  fseek(data_file, pstr->size*pstr->slice.num*sizeof(REAL), SEEK_CUR);
      }
      if ( DumpStream ) {
        v_tag = pvar->label;
	if (if_debug) printf ("Write Var: <%s> indx= %d\n", v_tag, j+1);
	rdwr_var (idf, j+1, bsize, v_tag, pvar, pbuf);
        sprintf (buf, "T\0", pcnt);
        odb_wrgr_(&idf, buf , pstr->time.data , strlen(buf));
      }
      else {
        v_tag = strcat(strcpy(v_tag, s_tag), pvar->label);
	if (if_debug) printf ("Write Var: <%s> indx= %d\n", v_tag, j+1);
	rdwr_var (idf, j+1, bsize, v_tag, pvar, pbuf);
      }}
    }
  }
  free (pbuf);
}


int do_stream (idf, pstr)
int idf;
STREAM *pstr;
{

  static  char s_cnt[4], sb[1000], vb[1000], tb[1000];
  static  int cnt;
  char *p1;

  /* make a stream number as 2-digit 0-padded number: */
  cnt++;
  sprintf(s_cnt, "%d\0", 100+cnt);

  /* make a stream base name: "TAG." */
  /*strcat(strcpy(sb, pstr->label), "."); changed -NHN*/
    strcat(strcpy(sb, pstr->label), "_");;

  if       (pstr->type == TP_STREAM) 
    def_stream  (idf, pstr, &s_cnt[1], sb, vb, tb);
  else 
    printf ("Unknown stream type! Stop!\n"), exit(-1);;

  write_stream (idf, pstr, &s_cnt[1], sb, vb); 
}

int do_grids (idf, cnt, nx, xx, ny, yy, nz, zz, nt, tt)
int    idf, nx, ny, nz, nt;
char  *cnt;
float *xx, *yy, *zz, *tt;
{
  char buf[16];
  int   i, itime, ind;
  REAL *tmp;
 
  tmp = (REAL *) malloc (nt * sizeof (REAL));

    if ( DumpVar || DumpStream) 
      sprintf (buf, "X\0", cnt);
    else
      sprintf (buf, "X_%s\0", cnt);
  def_grid(idf, buf, "Longitude", "degree_east", nx, xx);
  
    if ( DumpVar || DumpStream) 
      sprintf (buf, "Y\0", cnt);
    else
      sprintf (buf, "Y_%s\0", cnt);
  def_grid(idf, buf, "Latitude", "degree_north", ny, yy);
  
    if ( DumpVar || DumpStream) {
      sprintf (buf, "L\0", cnt);
      def_grid(idf, buf, "Level", "meters", nz, zz);}
    else {
      sprintf (buf, "L_%s\0", cnt);
      def_grid(idf, buf, "Level", "level", nz, zz);}
  ; 

    if ( DumpVar || DumpStream) 
      sprintf (buf, "T\0", cnt);
    else
      sprintf (buf, "T_%s\0", cnt);

      if (if_netcdf4grads) {
      for (i = 0; i <  nt; i++) tmp[i] = 365*tt[i]/12;
  def_grid(idf, buf, "Time", "days since 1960-01-01", nt, tmp);}
    else {
       if ( DumpVar || DumpStream ) 
  def_grid2(idf, buf, "Time", "months since 1960-1-1", nt, tt);
    else
  def_grid(idf, buf, "Time", "months since 1960-1-1", nt, tt);}
}


/* this is sort of nice but needs some extra work by Benno: 
 #define SET_LAT_LON_SCALES
*/

int def_grid (idf, name, lname, units, n, data)
int    idf, n;
char  *name, *lname, *units; 
float *data;
{
  int idv, sta, cnt;

    odb_dfgr_(&idf, name, &n, strlen(name));
    odb_setcattr_(&idf, name, "long_name", lname, strlen(name),9,strlen(lname));
    odb_setcattr_(&idf, name, "units", units,    strlen(name),5,strlen(units));
    odb_wrgr_(&idf, name, data, strlen(name)); 
}

int def_grid2 (idf, name, lname, units, n, data)
int    idf, n;
char  *name, *lname, *units; 
float *data;
{
  int idv, sta, cnt;

    odb_dftm_(&idf,name);
    odb_setcattr_(&idf, name, "long_name", lname, strlen(name),9,strlen(lname));
    odb_setcattr_(&idf, name, "units", units, strlen(name),5,strlen(units));
/*    odb_wrgr_(&idf, name, data, strlen(name)); */
}

int def_var (idf, name, lname, comp, cnt, xyz)
int    idf, comp;
char  *name, *lname, *xyz, *cnt;
{
  int idv, i, ndim = 4;
  char idim[4][8];

    if ( DumpVar || DumpStream) 
       for (i = 0; i < ndim; i++) sprintf (idim[i], "%c\0", xyz[i]);
    else
       for (i = 0; i < ndim; i++) sprintf (idim[i], "%c_%s\0", xyz[i], cnt);

    odb_dfvar_(&idf , &ndim, idim, name, 8, strlen(name));
    odb_setrattr_(&idf, name, "missing_value", &fmiss, strlen(name), 13);
   odb_setcattr_(&idf, name,"long_name",lname,strlen(name),9,strlen(lname));
}

int rdwr_var (idf, it, n, name, pvar, data)
int    idf, it, n;
char  *name;
VAR   *pvar;
float *data;
{
  static int sta[4], cnt[4];

  loc_rd (sizeof(REAL), n, data, data_file);

  odb_wr3v4_(&idf, &it, name, data, strlen(name));
}


int main (argc, argv)
int argc;
char **argv;
{
  static int ch, if_work=1;
  char space[200];
  int value;
  extern char *optarg;
  extern int   optind;

  nbegin = 1;
  nend = 999999;
  njump = 1;
  if (argc == 1) {
    printf ("%s: TIOS -> NetCDF conversion program\n\n", argv[0]);
    printf ("usage: %s [-tThGp][-f <data_file>][-o <out_file>][-s|v name]\n", argv[0]);
    exit (0);
  }
  else
    if_netcdf = 1;
    if_netcdf4grads = 0;
    while ((ch = getopt(argc, argv, "CGNThgtpf:v:V:o:s:")) != -1)
    switch (ch) {
    case 'h':
      printf ("%s: TIOS -> NetCDF conversion program\n\n", argv[0]);
      printf ("usage: %s [-thGp][-f <data_file>][-o <out_file>][-s|v name]\n", argv[0]);
      printf("\n\t-t - prints (only) the table of contents for <data_file>\n");
      printf("\t-T - prints (only) the names of streams in <data_file>\n");
      printf("\t-f <data_file> - Input file\n");
      printf("\t-p             - pick begin:end:jump for subsampling in time\n");
      printf("\t-o <out_file>  - Output file\n");
      printf("\t-s STREAM      - Dump this stream only\n");
      printf("\t-v VAR         - Dump this variable only\n");
      printf("\t-V VAR         - Dump this variable only, with full name\n");
      printf("\t-G             - Dump in GrADS-limited cdf format \n");
      printf("\t-g             - extra debugging info\n\n");
      printf("\t-h             - prints this help page\n\n");
      printf("defaults: <data_file>   : <%s>\n",   F_DATA);
      printf("          <out_file>   : <%s>\n\n", F_OUT);
      printf("Questions:\tnaomi@ldeo.columbia.edu\n");
      exit (0);
    case 'C':
      if_netcdf = 0;
      break;
    case 'G':
      if_netcdf4grads = 1;
      break;
    case 'O':
      if_old = 1;
      break;

    case 'v':
      DumpVar = strdup(optarg);
      break;
    case 'V':
      DumpFullVar = strdup(optarg);
      break;
    case 's':
      DumpStream = strdup(optarg);
      break;
    case 'p':
      Pick = 1;
      break;
    case 't': if_table = 1; if_work = 0;     break;
    case 'T': if_list = 1; if_work = 0;     break;
    case 'g': if_debug = 1;                  break;
    case 'f':
      strncpy(space, optarg, strlen(optarg));
      f_name = strdup(optarg);
      strcpy(space+strlen(optarg), ".indx");
      f_indx = strdup(space);
      if (!(indx_file = fopen((f_indx), "r"))) {
	printf ("cannot open data file: <%s>\n", f_indx);
	exit (-1);
      }
      strcpy(space+strlen(optarg), ".data");
      f_data = strdup(space);
      if (!(data_file = fopen((f_data), "r"))) {
	printf ("cannot open data file: <%s>\n", f_data);
	exit (-1);
      }
      if ( ! f_out) {
	strcpy(space+strlen(optarg), (if_netcdf ? ".nc" : ".cuf"));
	f_out = strdup(space);
      }
      break;

    case 'o': f_out = strdup(optarg); break;
    case '?':
      printf ("%s: TIOS -> NetCDF conversion program\n\n", argv[0]);
      printf ("usage: %s [-th][-f <data_file>][-o <out_file>][-s|v name]\n", argv[0]);
      exit (0);
      break;
    }

  if (!indx_file && !(indx_file = fopen((f_indx=F_INDX), "r"))) {
    printf ("cannot open data file: <%s>\n", F_INDX);
    exit (-1);
  }
  if (!data_file && !(data_file = fopen((f_data=F_DATA), "r"))) {
    printf ("cannot open data file: <%s>\n", F_DATA);
    exit (-1);
  }

  if ( ! f_out) f_out = F_OUT;

/*  if ( Pick ) {
     (void) printf("enter integer values: begin end step values   ");
     (void) fgets(space, sizeof(space), stdin);
     (void) sscanf(space,"%d %d %d", &nbegin, &nend, &njump);
  } */

  do_table ();
  if (if_work == 1) do_cdf();
}
/**************************************************************************/


