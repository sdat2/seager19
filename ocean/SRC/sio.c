/*********************************************************************
  This is  major text of TIOS library
  $Id: sio.c,v 1.1 1997/03/21 20:43:46 senya Exp senya $
*********************************************************************/
#include <stdio.h>
#include <string.h>
#include <values.h>
#include <ctype.h>
#include <stdlib.h>

#include <sys/types.h>
#include <malloc.h>

#include "tios.h"

static int sz_REAL, sz_UINT;

int idvar_tios_(a1,a2,a3) void *a1, *a2, *a3; {return tios_idvar_(a1,a2,a3);}

void loc_wr(size, num, buff, file)
int   size, num;
void *buff;
FILE *file;
{
  if (fwrite (buff, (size_t)size, (size_t)num, file) != num) {
    fprintf (stderr, "tios: disk-write error. Stop!\n");
    exit (-1);
  }
}
void loc_rd (size, num, buff, file)
int size, num;
void *buff;
FILE *file;
{
  if (fread (buff, (size_t)size, (size_t)num, file) != num) {
    fprintf (stderr, "tios: disk-read error. Stop!\n");
    exit (-1);
  }
}

/*.......................................................................*/
static int ipick (x0, n, xx)
/*.......................................................................*/
int n;
REAL x0, *xx;
{
  int i;

  if      (x0 <= xx[0]) 
    return 0;
  else if (x0 >= xx[n-1])
    return (n-1);
  else
    for (i = 1; i < n; i++) {
      if (x0 < xx[i]) 
	return (ABS(xx[i]-x0)<ABS(x0-xx[i-1])? i : (i-1));
    }
  return -1;
}
/*.......................................................................*/
static int ipick0 (x0, n, xx)
/*.......................................................................*/
int n;
REAL x0, *xx;
{
  int i;

  if      (x0 < xx[0]) 
    return -1;
  else if (x0 > xx[n-1])
    return -1;
  else
    for (i = 1; i < n; i++) {
      if (x0 < xx[i]) 
	return (ABS(xx[i]-x0)<ABS(x0-xx[i-1])? i : (i-1));
    }
  return n;
}
/*.......................................................................*/
static int read_array (arr)
/*.......................................................................*/
REAL **arr;
{
  REAL val[MAX_SLICES], tmp;
  int i, count = 0;
  char *word, *ch;

  if (!(word = strtok (NULL, " \t[")) || *word == ']') return 0; 

  do {
    if (!strncmp (word, "all", strlen("all"))) 
      return -1;
    else
      if (ch = strchr(word, ']')) {
	*ch = ' ';
	if (sscanf (word, "%f", &tmp) == 1) val[count++] = tmp;
	break;
      }
      else
	if (sscanf (word, "%f", &tmp) == 1) val[count++] = tmp;
  }
  while (word = strtok (NULL, " \t"));
    
  *arr = (REAL *) malloc(sz_REAL*count);
  for (i = 0; i < count; i++) *(*arr+i) = val[i];

  return count;
}

/*.......................................................................*/
static void pick_slice (pbuf, pstr, base1, n2, n3, map2, map3, co2, co3, pdat)
/*.......................................................................*/
int base1, n2, co2, *map2,
           n3, co3, *map3;
REAL *pbuf, *pdat;
STREAM *pstr;
{
  int j, k1, k2, k3, base;

  if (map3) {
    if (map2) {
      for (k1 = 0; k1 < pstr->slice.num; k1++) {
	base = *(pstr->slice.adr + k1) * base1;
	
	for (k2 = 0; k2 < n2; k2++) {
	  j = base + co2*map2[k2];
	  for (k3 = 0; k3 < n3; k3++) *pbuf++ = pdat[j + co3*map3[k3]];
	}
      }
    }
    else {
      for (k1 = 0; k1 < pstr->slice.num; k1++) {
	base = *(pstr->slice.adr + k1) * base1;
	
	for (k2 = 0; k2 < n2; k2++) {
	  for (k3 = 0; k3 < n3; k3++) *pbuf++ = pdat[base + co3*map3[k3]];
	  base += co2;
	}
      }
    }
  }
  else if (map2) {
    for (k1 = 0; k1 < pstr->slice.num; k1++) {
      base = *(pstr->slice.adr + k1) * base1;
      
      for (k2 = 0; k2 < n2; k2++) {
	j = base + co2*map2[k2];
	for (k3 = 0; k3 < co3*n3; k3+=co3) *pbuf++ = pdat[j + k3];
      }
    }
  }
  else {
    for (k1 = 0; k1 < pstr->slice.num; k1++) {
      base = *(pstr->slice.adr + k1) * base1;
      
      for (k2 = 0; k2 < n2; k2++) {
	for (k3 = 0; k3 < co3*n3; k3+=co3) *pbuf++ = pdat[base + k3];
	base += co2;
      }
    }
  }
}
/*.......................................................................*/
static void pick_sprite (pbuf,pstr,base1,n2,n3,map2,map3,co2,co3,pdat)
/*.......................................................................*/
int base1, n2, co2, *map2,
           n3, co3, *map3;
REAL *pbuf, *pdat;
STREAM *pstr;
{
  int j, m, r, k1, k2, k3, base, NXY, NMAP, *cmp;
  REAL flag;
  VAR  *pvar = pstr->vars.adr[0];
  
  flag = pvar->flag;
  cmp  = pvar->map->cmp;
  NMAP = pvar->map->NMAP; 
  NXY  = pvar->grid->NX * pvar->grid->NY;

  if (map3) {
    if (map2) {
      for (k1 = 0; k1 < pstr->slice.num; k1++) {
	base = *(pstr->slice.adr + k1) * base1;
	
	for (k2 = 0; k2 < n2; k2++) {
	  j = base + co2*map2[k2];
	  for (k3 = 0; k3 < n3; k3++) {
	    m = j + co3*map3[k3];
	      *pbuf++ = (r = cmp[m % NXY]) ? pdat[r + (m/NXY)*NMAP - 1] : flag;
	  }
	}
      }
    }
    else {
      for (k1 = 0; k1 < pstr->slice.num; k1++) {
	base = *(pstr->slice.adr + k1) * base1;
	
	for (k2 = 0; k2 < n2; k2++) {
	  for (k3 = 0; k3 < n3; k3++) {
	    m = base + co3*map3[k3];
	    *pbuf++ = (r = cmp[m % NXY]) ? pdat[r + (m/NXY)*NMAP - 1] : flag;
	  }	      
	  base += co2;
	}
      }
    }
  }
  else if (map2) {
    for (k1 = 0; k1 < pstr->slice.num; k1++) {
      base = *(pstr->slice.adr + k1) * base1;
      
      for (k2 = 0; k2 < n2; k2++) {
	j = base + co2*map2[k2];
	for (k3 = 0; k3 < co3*n3; k3+=co3) {
	  m = j + k3; 
	  *pbuf++ = (r = cmp[m % NXY]) ? pdat[r + (m/NXY)*NMAP - 1] : flag;
	}
      }
    }
  }
  else {
    for (k1 = 0; k1 < pstr->slice.num; k1++) {
      base = *(pstr->slice.adr + k1) * base1;
      
      for (k2 = 0; k2 < n2; k2++) {
	for (k3 = 0; k3 < co3*n3; k3+=co3) {
	  m = base + k3;
	  *pbuf++ = (r = cmp[m % NXY]) ? pdat[r + (m/NXY)*NMAP - 1] : flag;
	}
	base += co2;
      }
    }
  }
}
/*.......................................................................*/
static void fill_buffer (pstr, pbuf, pdat)
/*.......................................................................*/
REAL   *pbuf, *pdat;
STREAM *pstr;
{
  int NX, NY, NXY, MX, MY, MZ, *mapx, *mapy, *mapz;
  GR_MAP *pgm  = pstr->gmap;
  VAR    *pvar = pstr->vars.adr[0];

  mapx = mapy = mapz = NULL;

  MX = NX = pvar->grid->NX; 
  MY = NY = pvar->grid->NY; 
  MZ =      pvar->grid->NZ;
  NXY = NX * NY;

  if (pgm) {
    if (pgm->MX) MX = pgm->MX, mapx = pgm->mapx;
    if (pgm->MY) MY = pgm->MY, mapy = pgm->mapy;
    if (pgm->MZ) MZ = pgm->MZ, mapz = pgm->mapz;
  }

  switch (pstr->slice.id) {
  case 'X': case 'x':
    if (pvar->map)
      pick_sprite (pbuf,pstr, 1,MZ,MY, mapz, mapy, NXY,NX, pdat);
    else
      pick_slice (pbuf, pstr, 1,  MZ, MY, mapz, mapy, NXY,NX, pdat);
    break;
    
  case 'Y': case 'y':
    if (pvar->map)
      pick_sprite (pbuf,pstr, NX, MZ,MX, mapz, mapx, NXY,1,pdat);
    else
      pick_slice (pbuf, pstr, NX,  MZ, MX, mapz, mapx, NXY,1, pdat);
    break;
    
  case 'Z': case 'z':
    if (pvar->map)
      pick_sprite (pbuf,pstr, NXY, MY,MX, mapy, mapx, NX,1,pdat);
    else
      pick_slice (pbuf, pstr, NXY, MY, MX, mapy, mapx, NX, 1, pdat);
    break;
  }
}
/*.......................................................................*/
static void fill_cluster (ps, pbuf, pdat)
/*
 * fills the cluster buffer with data. 
 * Works both for compressed and uncompressed variables.
 *.......................................................................*/
REAL   *pbuf, *pdat;
STREAM *ps;
{
  int  i, j, k, nchunk, base;
  REAL flag;
  VAR *pvar = ps->vars.adr[0];
  
  flag = pvar->flag;

  if ( pvar->map ) {
    nchunk = pvar->map->NMAP;
    for (k = 0; k < ps->slice.num; k++) {

      base = ps->slice.adr[k] * nchunk;
      for (i = 0; i < ps->cluster.n; i++) {
	*pbuf++ = (j = ps->cluster.indx[i]) ? pdat[base + j] : flag ; 
	if (j >= nchunk) printf ("i=%d j=%d !!!\n", i, j);
      }
    }
  }
  else {
    nchunk = pvar->grid->NX * pvar->grid->NY;
    for (k = 0; k < ps->slice.num; k++) {

      base = ps->slice.adr[k] * nchunk;
      for (i = 0; i < ps->cluster.n; i++) 
	*pbuf++ = pdat[base + ps->cluster.indx[i]]; 
    }
  }
}

/*.......................................................................*/
static void write_stream_data (pstr, time)
/*.......................................................................*/
STREAM *pstr; 
REAL   *time;
{
  static REAL *buff;
  int buf_size, i;

  buf_size = (pstr->type == TP_STREAM ?	pstr->size : pstr->cluster.n);
  buf_size *= pstr->slice.num;
  if (! buf_size) printf("tios: zero buffer size! Stop.\n"), exit(-1); 

  buff = (REAL *) calloc(buf_size, sizeof(REAL) );

  pstr->time.count++;
  pstr->curr->time = *time;
  pstr->curr->val = ftell (data_file);

  for (i = 0; i < pstr->vars.count; i++) {
    if (pstr->type == TP_STREAM) 
      fill_buffer (pstr, buff, pstr->vars.adr[i]->var);
    else
      fill_cluster (pstr, buff, pstr->vars.adr[i]->var);

    loc_wr (sz_REAL, buf_size, (void *)buff, data_file);
  }
  /* allocate a new time-chain element; make it current: */
  pstr->curr = (pstr->curr->next = (ADDR *) calloc(1, sizeof(struct ADDR)));

  fflush (data_file);
  free (buff);
}

/*.......................................................................*/
static void wr_as_real (n, idat, file)
/*.......................................................................*/
int n, *idat;
FILE *file;
{
  int i; 
  REAL *real;
  real = (REAL *)malloc (sz_REAL*n);
  for (i = 0; i < n; i++) real[i] = (REAL)idat[i];
  loc_wr (sz_REAL, n, (void *)real, file);
  free (real);
}

/*.......................................................................*/
static void lwrite_tios (file)
/*.......................................................................*/
FILE *file;
{
  static REAL tmp[NEW_TIOS_HEADER];

  if ( tios.version ) {
    loc_wr (1, 8, Magic, file);
    tmp[0] = (REAL)tios.version;
    loc_wr (sz_REAL, 1,  (void *)tmp, file);
  }

  tmp[0]  = (REAL)tios.grid_count; 
  tmp[1]  = (REAL)tios.gmap_count;
  tmp[2]  = (REAL)tios.var_count;
  tmp[3]  = (REAL)tios.str_count;
  tmp[4]  = (REAL)tios.time_begin; 
  tmp[5]  = (REAL)tios.time_end;
  tmp[6]  = (REAL)tios.addr_start; 
  tmp[7]  = (REAL)tios.addr_grids;
  tmp[8]  = (REAL)tios.addr_vars;  
  tmp[9]  = (REAL)tios.addr_strs;  
  tmp[10] = (REAL)tios.addr_end;

  tmp[11] = (REAL)tios.clust_count;
  tmp[12] = (REAL)tios.pole_alpha;
  tmp[13] = (REAL)tios.pole_beta;
  tmp[14] = (REAL)tios.pole_gamma;

  if ( tios.version == 0) 
    loc_wr (sz_REAL, OLD_TIOS_HEADER,  (void *)tmp, file);
  else
    loc_wr (sz_REAL, NEW_TIOS_HEADER,  (void *)tmp, file);

  loc_wr (1, TITLE_TIOS_LEN, (void *)tios.tios_name, file);
}
/*.......................................................................*/
static void restart_tios (file)
/*.......................................................................*/
FILE *file;
{
  static REAL tmp[NEW_TIOS_HEADER];
  static char magic[8];

  loc_rd (1, 8, magic, file);
  if ( ! strncmp(magic, Magic, 8) ) { /* new version */
    loc_rd (sz_REAL, 1, (void *)tmp, file);
    tios.version = (int)tmp[0];
  }
  else {                              /* an old 0-version */
    fseek (file, 0, SEEK_SET);
    tios.version = 0;
  }

  if ( tios.version == 0) 
    loc_rd (sz_REAL, OLD_TIOS_HEADER,  (void *)tmp, file);
  else
    loc_rd (sz_REAL, NEW_TIOS_HEADER,  (void *)tmp, file);

  tios.time_begin =      tmp[4];
  tios.addr_start = (int)tmp[6];
  tios.addr_grids = (int)tmp[7];
  tios.addr_vars  = (int)tmp[8];
  tios.addr_strs  = (int)tmp[9];
  tios.addr_end   = (int)tmp[10];

  loc_rd (1, TITLE_TIOS_LEN, tios.tios_name, file);
}
/*.......................................................................*/
static void lwrite_grid (pg)
/*.......................................................................*/
GRID *pg;
{
  static REAL tmp[3];

  tmp[0] = (REAL)(pg->NX);   
  tmp[1] = (REAL)(pg->NY);
  tmp[2] = (REAL)(pg->NZ);   
  loc_wr (sz_REAL, 3, (void *)tmp, indx_file);

  if (pg->NX) loc_wr(sz_REAL, pg->NX, (void *)pg->xx, indx_file);
  if (pg->NY) loc_wr(sz_REAL, pg->NY, (void *)pg->yy, indx_file);
  if (pg->NZ) loc_wr(sz_REAL, pg->NZ, (void *)pg->zz, indx_file);
}
/*.......................................................................*/
static void lwrite_gmap (pgm)
/*.......................................................................*/
GR_MAP *pgm;
{
  static REAL tmp[4];

  tmp[0] = (REAL)(pgm->MX); tmp[1] = (REAL)(pgm->MY);
  tmp[2] = (REAL)(pgm->MZ); tmp[3] = (REAL)(pgm->grid->grid_nu);  
  loc_wr (sz_REAL, 4, (void *)tmp, indx_file);

  if (pgm->MX) wr_as_real(pgm->MX, pgm->mapx, indx_file);
  if (pgm->MY) wr_as_real(pgm->MY, pgm->mapy, indx_file);
  if (pgm->MZ) wr_as_real(pgm->MZ, pgm->mapz, indx_file);
}
/*.......................................................................*/
static void lwrite_var (pv)
/*.......................................................................*/
VAR *pv;
{
  static REAL tmp[3];

  tmp[0] = (REAL)(pv->strs.count);  
  tmp[1] = (REAL)(pv->grid->grid_nu);
  tmp[2] = pv->flag;
  loc_wr (sz_REAL, 3, (void *)tmp, indx_file);

  if (pv->strs.count) wr_as_real (pv->strs.count, pv->strs.num, indx_file);
  loc_wr (1, STR_VAR_LEN, pv->label, indx_file);
}
/*.......................................................................*/
static void lwrite_stream (ps)
/*.......................................................................*/
STREAM *ps;
{
  int nsz;
  static REAL tmp[9];

  if ( tios.version ) {
    tmp[0] = (REAL)ps->type;
    loc_wr (sz_REAL, 1, (void *)tmp, indx_file);
  }

  tmp[0] = (REAL)(ps->time.count);  
  tmp[1] = ps->time.sta;  
  tmp[2] = ps->time.sta + ps->time.step * (REAL)(ps->time.count - 1);  
  tmp[3] = ps->time.step; 
  tmp[4] = (REAL)(ps->slice.id);
  tmp[5] = (REAL)(ps->slice.num); 
  tmp[6] = (REAL)(ps->gmap ? ps->gmap->gmap_nu : 0);
  tmp[7] = (REAL)(ps->size); 
  tmp[8] = (REAL)(ps->vars.count);

  loc_wr (sz_REAL, 9, (void *)tmp, indx_file);
  loc_wr (sz_REAL, ps->slice.num, (void *)ps->slice.val, indx_file);
  wr_as_real ((int)(ps->vars.count), ps->vars.num, indx_file); 

  loc_wr (1, LABEL_STREAM_LEN, ps->label, indx_file);

  if ( ps->type == TP_CLUSTER ) { /* for a cluster-type output grids: */

    nsz = ps->cluster.n;

    tmp[0] = (REAL)(nsz);
    tmp[1] = (REAL)(ps->cluster.nid);  
    tmp[2] = (REAL)(ps->cluster.ntag);  
    tmp[4] = (REAL)0;
    loc_wr (sz_REAL, 4, (void *)tmp, indx_file);
    loc_wr (sz_REAL, nsz, (void *)ps->cluster.px, indx_file);
    loc_wr (sz_REAL, nsz, (void *)ps->cluster.py, indx_file);
    if (ps->cluster.nid)
      loc_wr (sz_UINT, nsz*ps->cluster.nid, (void *)ps->cluster.pid, indx_file);
    if (ps->cluster.ntag)
      loc_wr (1, nsz * ps->cluster.ntag, (void *)ps->cluster.ptag, indx_file);
  }

  fflush(indx_file);
}
/*.......................................................................*/
static void restart_stream (pstr)
/*.......................................................................*/
STREAM *pstr;
{
  int i;
  ADDR *paddr;
  UINT *pint;
  REAL *pflt;
  static REAL tmp[9];

  if ( tios.version ) {
    loc_rd (sz_REAL, 1, (void *)tmp, indx_file);
    if ( pstr->type != (int)tmp[0] ) {
      printf ("TIOS: in restart_stream() old type(%d) != new type(%d). Stop!\n",
	      (int)tmp[0], pstr->type);
      exit (-1);
    }
  }

  loc_rd (sz_REAL, 9, (void *)tmp, indx_file);
  pstr->time.count = (int)tmp[0];
  pstr->time.sta   =      tmp[1];
  pstr->time.next  =      tmp[2];

  FSKIP (indx_file, (((int)tmp[5] + (int)tmp[8])*sz_REAL + LABEL_STREAM_LEN));

  if ( pstr->type == TP_CLUSTER ) {                /* for a cluster */
    loc_rd (sz_REAL, 4, (void *)tmp, indx_file);
    if ( pstr->cluster.n != (int)tmp[0] ) {
      printf ("TIOS: in restart_stream() old N(%d) != new N(%d). Stop!\n",
	      (int)tmp[0], pstr->cluster.n);
      exit (-1);
    }
    if ( pstr->cluster.nid != (int)tmp[1] || 
	pstr->cluster.ntag != (int)tmp[2] ) {
      printf ("TIOS: in restart_stream() CLUSTER re-defined. Stop!\n");
      exit (-1);
    }
    FSKIP (indx_file, pstr->cluster.n * sz_REAL);
    if (pstr->cluster.nid) FSKIP (indx_file, pstr->cluster.n * sz_UINT);
    if (pstr->cluster.ntag) FSKIP(indx_file,pstr->cluster.n*pstr->cluster.ntag);
  }

  paddr = pstr->addr;

  if (pstr->time.count) {
    pint = (UINT *)malloc(sz_UINT*pstr->time.count);
    pflt = (REAL *)malloc(sz_REAL*pstr->time.count);

    loc_rd (sz_UINT, pstr->time.count, (void *)pint, indx_file);
    loc_rd (sz_REAL, pstr->time.count, (void *)pflt, indx_file);

    for (i = 0; i < pstr->time.count; i++) {
      paddr->val  = pint[i];
      paddr->time = pflt[i];
      paddr = (paddr->next = (ADDR *) calloc(1, sizeof (struct ADDR)));
    }
    free(pint);
    free(pflt);
  }

  pstr->curr = paddr;
}
/*.......................................................................*/
static void do_gmap (MM, map, rang, NN, xx)
/*.......................................................................*/
RG_ST *rang;
int NN, *MM, **map;
REAL *xx;
{
  int i, st, count = 0, *pm;
  REAL first, last, delt;

  if (!rang->num || NN == 1)
    return;

  else { 
    first = (rang->sta == AUTO_FLAG) ? xx[0]    : rang->sta;
    last  = (rang->end == AUTO_FLAG) ? xx[NN-1] : rang->end;
    
    if (rang->num < 0) {
      st = ((rang->step  == AUTO_FLAG) ? 1 : (int)(rang->step));
      
      for (i = 0; i < NN; i += st) {
	if (xx[i] >= first && xx[i] <= last) count++;
      }
      
      *MM = count;
      *map = pm = (int *)malloc(count * sizeof(int));
      
      for (count = i = 0; i < NN; i += st) {
	if (xx[i] >= first && xx[i] <= last) pm[count++] = i;
      }
    }
    else {
      *MM = count = MIN(rang->num, NN);
      
      *map = pm = (int *)malloc(count * sizeof(int));
      delt = (last - first) / (REAL)(count - 1);
      
      for (i = 0; i < count; i++) {
	pm[i] = ipick (first, NN, xx);
	first += delt;
      }
    }
  }
}
/*.......................................................................*/
static void regrid_stream (nran, pstr, pvar)
/*.......................................................................*/
STREAM *pstr;
VAR    *pvar;
int     nran;
{
  int ix, iy, i, j, ntag;
  GRID   *pg = pvar->grid;
  GR_MAP *pgm;
  RANGE  *pr;

  if ( ! nran ) {
    pstr->time.sta  = pstr->time.end  = AUTO_FLAG;
    pstr->time.step = 1.;
    pstr->time.fmt  = 0;
    return;
  }
  else {
    pr = tios.ranges;
    while (pr && pr->range_nu != nran) pr = pr->next;

    if (!pr) {
      printf("TIOS: Not defined RANGE %d used in stream %d for <%s>\n", 
	     nran, pstr->stream_nu, pvar->label);
      exit (-1);
    }

    pstr->time.sta  = pr->t.sta;
    pstr->time.end  = pr->t.end;
    pstr->time.step = pr->t.step;

    pstr->time.fmt = pr->t.fmt;
    pstr->time.absolute = pr->t.absolute;

    if ( pr->type == TP_CLUSTER ) { 
/*......................................................................
 * for a cluster-type streams: allocates the memory for the cluster;    
 * copies the cluster grid information from the RANGE
 *.....................................................................*/
      pstr->type = TP_CLUSTER; 

      pstr->cluster.px   = (REAL *)calloc(pr->cluster.n, sizeof(REAL));
      pstr->cluster.py   = (REAL *)calloc(pr->cluster.n, sizeof(REAL));
      pstr->cluster.indx = (UINT *)calloc(pr->cluster.n, sizeof(UINT));
      if (pr->cluster.nid) 
	pstr->cluster.pid = (UINT *)calloc(pr->cluster.n, sizeof(UINT));
      if (ntag = pr->cluster.ntag) 
	pstr->cluster.ptag =(char *)calloc(1, pr->cluster.n * ntag);
      
      for (i = j = 0; i < pr->cluster.n; i++) {

	ix = ipick0(pr->cluster.px[i], pg->NX, pg->xx);
	if (ix == -1) { 
	  printf("TIOS: cluster point (%g,%g) is out of X-range [%g : %g]. Ignored...!\n", 
	  pr->cluster.px[i], pr->cluster.py[i],pg->xx[0],pg->xx[pg->NX-1]); 
	  continue; /* skip this point */
	}
	iy = ipick0(pr->cluster.py[i], pg->NY, pg->yy);
	if (iy == -1) { 
	  printf("TIOS: cluster point (%g,%g) is out of Y-range [%g : %g]. Ignored...!\n", 
	  pr->cluster.px[i], pr->cluster.py[i],pg->yy[0],pg->yy[pg->NY-1]); 
	  continue; /* skip this point */
	}
	pstr->cluster.px[j] = pg->xx[ix];
	pstr->cluster.py[j] = pg->yy[iy];
	pstr->cluster.indx[j] =
	  ( pvar->map ? pvar->map->cmp[iy * pg->NX + ix] : iy * pg->NX + ix );

	if (pr->cluster.nid) pstr->cluster.pid[j] = pr->cluster.pid[i];
	if ( ntag )
	  memcpy(&pstr->cluster.ptag[j*ntag], &pr->cluster.ptag[i*ntag], ntag);
	
	j++;
      }
      if (pstr->cluster.n = j) {  /* for a non-empty cluster: */
	pstr->cluster.type = pr->cluster.type;
	pstr->cluster.nid  = pr->cluster.nid;
	pstr->cluster.ntag = pr->cluster.ntag;
      }
      else {
	printf("TIOS: Empty Sream/Cluster <%s>. Stop.\n", pstr->label);
	exit (-1);
      }

      return;
    }
    
             /* if there are no sub-sampling in any of the grids - return */
    if (!pr->x.num && !pr->y.num && !pr->z.num) return;
    
    pgm = tios.gr_maps;
    while (pgm->next && (pgm->grid != pg || pgm->range != pr))
      pgm = pgm->next;
    
    if ( !pgm->next ) { /* a new grid-map should be added: */
      do_gmap (&(pgm->MX), &(pgm->mapx), &(pr->x), pg->NX, pg->xx);
      do_gmap (&(pgm->MY), &(pgm->mapy), &(pr->y), pg->NY, pg->yy);
      do_gmap (&(pgm->MZ), &(pgm->mapz), &(pr->z), pg->NZ, pg->zz);
      
      pgm->gmap_nu = ++tios.gmap_count;
      pgm->grid    = pg;
      pgm->range   = pr;
      pstr->gmap   = pgm;
      pgm->next    = (GR_MAP *)calloc(1, sizeof(struct GR_MAP)); 
    }
    else                /* found matching previously defined grid-map */
      pstr->gmap = pgm;
  }
}

/*.......................................................................*/
static void set_slice (pstr, pvar)
/*.......................................................................*/
STREAM *pstr;
VAR    *pvar;
{
  REAL *xx;
  int i, nn;
  
  switch (pstr->slice.id) {
  case 'X':  case 'x':
    nn = pvar->grid->NX;
    xx = pvar->grid->xx;
    pstr->size = 
      ((pstr->gmap && pstr->gmap->MY)? pstr->gmap->MY : pvar->grid->NY) *
      ((pstr->gmap && pstr->gmap->MZ)? pstr->gmap->MZ : pvar->grid->NZ);
    break;
  case 'Y':  case 'y':
    nn = pvar->grid->NY;
    xx = pvar->grid->yy;
    pstr->size = 
      ((pstr->gmap && pstr->gmap->MX)? pstr->gmap->MX : pvar->grid->NX) *
      ((pstr->gmap && pstr->gmap->MZ)? pstr->gmap->MZ : pvar->grid->NZ);
    break;
  case 'Z':  case 'z':
    nn = pvar->grid->NZ;
    xx = pvar->grid->zz;
    pstr->size = 
      ((pstr->gmap && pstr->gmap->MX)? pstr->gmap->MX : pvar->grid->NX) *
      ((pstr->gmap && pstr->gmap->MY)? pstr->gmap->MY : pvar->grid->NY);
    break;
  default:
    printf ("TIOS: wrong axe name <%c> has been specified\n", pstr->slice.id);
    exit (-1);
  }

  if (pstr->slice.num == -1) {
    pstr->slice.num = nn;
    pstr->slice.val = xx;
  }

  pstr->slice.adr = (int *) malloc(pstr->slice.num * sizeof(int));
  for (i = 0; i < pstr->slice.num; i++) /* a black magic loop to puzzle you: */
    pstr->slice.val[i] = 
      xx[pstr->slice.adr[i] = 
	 ipick(pstr->slice.val[i], nn, xx)];
}
/*.......................................................................*/
static void set_define ()
/*.......................................................................*/
{
  char *word;

  if (!(word = strtok (NULL, " \t\n")))
    {printf ("TIOS warning: empty DEFINE\n"); return;}
  
  if (!strncmp (word, "TIME_FMT", strlen("TIME_FMT"))) {
    word = strtok (NULL, " \t\n[]");
    switch  (*word) {
    case 's':
      tios.time_fmt = 0;
      break;
    case 'd':
      tios.time_fmt = 1;
      break;
    case 'm':
      tios.time_fmt = 2;
      break;
    case 'y':
      tios.time_fmt = 3;
      break;
    default:
      tios.time_fmt = 2;
      break;
    }
  }
  else if (!strncmp (word, "LABEL", strlen("LABEL"))) 
    strcpy(tios.tios_name,
	   strtok(strpbrk(word+strlen("LABEL")+1,"[")+1, "]\n") );

  else if (!strncmp (word, "DUMP_OUTPUT", strlen("DUMP_OUTPUT"))) 
    tios.cntrl.everystep = 1;

  else if (!strncmp (word, "DEBUG", strlen("DEBUG"))) 
    tios.cntrl.debug = 1;
}
/*------------------------------------------------------------------*/
static void read_range(rang, word)
/*------------------------------------------------------------------*/
RG_ST *rang;
char  *word;
{
  int i=0;
  char fmt, grid_name = *word;

  while ( word = strtok (NULL, " \t:;" ) ) {
    i++;

    if (*word == '*' ||
	!strncmp(word, (i-2)?FRST_WORD:LAST_WORD, strlen(word)) ) continue;

    if (*word == ']') break;

    switch (i) {
      case 1: sscanf (word, "%f", &(rang->sta)); rang->num = -1; break;
      case 2: sscanf (word, "%f", &(rang->end)); rang->num = -1; break;
      case 3: 
	if (*word == '(') sscanf (word+1, "%d", &(rang->num)); 
	else              sscanf (word,   "%f", &(rang->step)); 
	break;
      case 4: 
	if (!strncmp(word, "fmt=", strlen("fmt="))) {
	  sscanf (&word[4], "%c", &fmt);
	  switch  (fmt) {
	    case 's': case 'S': rang->fmt = 0; break;
	    case 'd': case 'D': rang->fmt = 1; break;
	    case 'm': case 'M': rang->fmt = 2; break;
	    case 'y': case 'Y': rang->fmt = 3; break;
	    default:  rang->fmt = tios.time_fmt; break;
	  }
	}
      case 5:
	if (!strncmp(word, "abs", strlen("abs")) ) rang->absolute = 1;      
      default: break;
    }

    if ( strchr(word, ']') ) break;    
  }

  if (tios.cntrl.debug) {
    printf ("[%c:", grid_name);
    (rang->sta == (REAL)AUTO_FLAG) ? printf ("*:"): printf ("%g:",rang->sta); 
    (rang->end == (REAL)AUTO_FLAG) ? printf ("*:"): printf ("%g:",rang->end);
    (rang->step == (REAL)AUTO_FLAG) ? printf ("*"): printf ("%g",rang->step);
    switch (rang->fmt) {
     case 0:  fmt = 's'; break; 
     case 1:  fmt = 'd'; break; 
     case 2:  fmt = 'm'; break; 
     case 3:  fmt = 'y'; break; 
     default: fmt = '?'; break; 
    }
    if (grid_name == 'T') printf ("; fmt=%c; abs=%d", fmt, rang->absolute);
    printf ("]");
  }
}

/*.......................................................................*/
static void set_range (line0, word)
/*.......................................................................*/
char *line0, *word;
{
  char *ic1 = line0;
  static char *line;
  static int lenline;
  
  if (!(word = strtok (NULL, " \t\n")))
    {printf ("TIOS: Missed RANGE ID\n");exit (-1);}
  sscanf (word, "%d", &(rang_curr->range_nu));

  rang_curr->type = TP_STREAM;

  rang_curr->x.sta = rang_curr->x.end = rang_curr->x.step = 
  rang_curr->y.sta = rang_curr->y.end = rang_curr->y.step = 
  rang_curr->z.sta = rang_curr->z.end = rang_curr->z.step = 
  rang_curr->t.sta = rang_curr->t.end = (REAL)AUTO_FLAG;
  rang_curr->t.step = 1.; 
  rang_curr->t.fmt  = tios.time_fmt;
  rang_curr->t.absolute = 0;

  if (tios.cntrl.debug) printf ("Debug: Range #%1d :", rang_curr->range_nu);

  while ( ic1 = strchr(ic1, '[') ) {

    if (lenline < strlen(ic1)) lenline = strlen(line = strdup(ic1));
    else                                 strcpy(line, ic1);

    word = strtok (line, " \t\n[");

    if ( !(ic1 = strchr(ic1, ']')) )
      {printf ("TIOS: \"]\" expected for <%c> RANGE\n", *word);exit (-1);}

    switch (*word) {

    case 'X': case 'x':
      read_range (&(rang_curr->x), word);
      break;
    case 'Y': case 'y':
      read_range (&(rang_curr->y), word);
      break;
    case 'Z': case 'z':
    case 'L': case 'l':
      read_range (&(rang_curr->z), word);
      break;
    case 'T': case 't':
      read_range (&(rang_curr->t), word);
      break;

    default:
      break;
    }
  }
  if (tios.cntrl.debug) printf ("\n");

  rang_curr->range_nu = ++tios.rang_count;
  rang_curr = (rang_curr->next = (RANGE *)calloc (1, sizeof(struct RANGE)));
}

/*.......................................................................*/
static void set_cluster (line0)
/*.......................................................................*/
char *line0;
{
  char   *word, *ic1 = line0, *pc1, *pc2, *fstation = (char *)0;
  static char *line;
  static int lenline;
  struct CLUSTER *pc;
  FILE *file;
  int i;
  
/*.......................................................................
 * parse the line : CLUSTER number size [T .....] [file= ] [id=] [tags=]
 *......................................................................*/

  if (!(word = strtok (NULL, " \t[]\n")))
    {printf ("TIOS: Missed CLUSTER ID\n");exit (-1);}
  sscanf (word, "%d", &(rang_curr->range_nu));

  if (!(word = strtok (NULL, " \t[]\n")))
    {printf ("TIOS: Missed CLUSTER size\n");exit (-1);}
  sscanf (word, "%d", &(rang_curr->cluster.n));

  if (rang_curr->cluster.n <= 0) 
    {printf ("TIOS: invalid CLUSTER size: %d\n",rang_curr->cluster.n);exit(-1);}

  rang_curr->type = TP_CLUSTER;
  pc = &rang_curr->cluster;

  rang_curr->x.sta = rang_curr->x.end = rang_curr->x.step = 
  rang_curr->y.sta = rang_curr->y.end = rang_curr->y.step = 
  rang_curr->z.sta = rang_curr->z.end = rang_curr->z.step = 
  rang_curr->t.sta = rang_curr->t.end = (REAL)AUTO_FLAG;
  rang_curr->t.step = 1.; 
  rang_curr->t.fmt  = tios.time_fmt;
  rang_curr->t.absolute = 0;

  if (tios.cntrl.debug) printf ("Debug: Cluster #%1d :", rang_curr->range_nu);

  while ( ic1 = strchr(ic1, '[') ) { 

    if (lenline < strlen(ic1)) lenline = strlen(line = strdup(ic1));
    else                                 strcpy(line, ic1);

    word = strtok (line, " \t\n[=");

    if ( !(ic1 = strchr(ic1, ']')) )
      {printf ("\nTIOS: \"]\" expected for CLUSTER <%s>\n", line0);exit (-1);}

    if      ( !strncmp(word, "T", 1) || !strncmp(word, "t", 1) ) 
      read_range (&(rang_curr->t), word);
    else {
      printf ("\nTIOS: invalid CLUSTER parameter! Stop.\n--> %s\n", line0);
      exit (-1);
    }
  }

  while ( word = strtok (NULL, " \t\n=")) {/* parsing the filename,IDs,Tags: */
    if      ( !strncasecmp(word, "file", 4) && (word = strtok (NULL, " \t\n")))
      fstation = strdup(word);
    else if ( !strncasecmp(word, "id", 2) )
      pc->nid = 1;
    else if ( !strncasecmp(word, "tag", 3) && (word = strtok (NULL, " \t\n")))
      pc->ntag = atoi(word);
  } 

  pc->px = (REAL *)calloc(pc->n, sizeof(REAL));
  pc->py = (REAL *)calloc(pc->n, sizeof(REAL));
  if (pc->nid)  pc->pid  = (UINT *)calloc(pc->n, sizeof(UINT));
  if (pc->ntag) pc->ptag = (char *)calloc(pc->n * pc->ntag, sizeof(char));

  if ( fstation ) {
    if ( !(file = fopen(fstation, "r")) ) {
      printf ("\nTIOS: cannot open station file <%s>! Stop.\n", fstation);
      exit(-1);
    }
    else { 
      fseek(file,0,SEEK_SET);
      if (tios.cntrl.debug) printf(" reading file <%s>:", fstation);
    }
  }
  else file = tios_file;

  i = 0;
  while ( i < pc->n ) {                    /* parsing the stations file */
    if ( ! fgets(line, TITLE_TIOS_LEN, file) ){
      printf ("\nTIOS: not enough data for CLUSTER <%d>\n",rang_curr->range_nu);
      exit(-1);
    }

    if (*line == '%') continue;           /* skip the comment */
    
    if (sscanf (line, "%f %f", &(pc->px[i]), &(pc->py[i])) != 2) {
      printf ("\nTIOS: incomplete statiion data! Stop.\n%d--> %s\n",i,line);
      exit(-1);
    }
    if (pc->nid && (sscanf (line, "%*f %*f %d", &(pc->pid[i])) != 1) ) {
      printf ("\nTIOS: incomplete statiion ID data! Stop.\n%d--> %s\n",i,line);
      exit(-1);
    }
    if (pc->ntag) {
      if ( !(pc1 = strpbrk(line, "\'\"([")) || !(pc2 = strrchr(line, *pc1))) {
     printf ("\nTIOS: incomplete statiion NAME data! Stop.\n%d--> %s\n",i,line);
	exit(-1);
      }
      *pc2 = '\0'; 
      strncpy (&(pc->ptag[i*pc->ntag]), pc1+1, pc->ntag); 
    }
    i++;
  }
  if ( fstation ) fclose(file);

  if (tios.cntrl.debug) {
    printf ("\n");
    for (i = 0; i < pc->n; i++) {
      printf ("station: %3d x=%g y=%g", i+1, pc->px[i],pc->py[i]);
      if (pc->nid)  printf ("\tID=%d", pc->pid[i]);
      if (pc->ntag) printf ("\tname=\"%s\"", &pc->ptag[i*pc->ntag]);
      printf ("\n");
    }
  }

  rang_curr->range_nu = ++tios.rang_count;
  rang_curr = (rang_curr->next = (RANGE *)calloc (1, sizeof(struct RANGE)));
}

/*.......................................................................*/
static VAR *find_var_by_name (name)
/*.......................................................................*/
char *name;
{
  VAR *pvar = tios.vars;
  
  while (pvar->next && strcmp(name, pvar->label))
    pvar = pvar->next; 

  return (pvar->next ? pvar : pvar->next);
}
/*.......................................................................*/
static void set_stream(line)
/*.......................................................................*/
char *line;     
{
  VAR    *pvar;
  REAL   *vslice;
  int    i, nran, nslice;
  char   *word, label[LABEL_STREAM_LEN], islice;

  if (!(word = strtok (NULL, " \t\n"))) return;   /* get LABEL */
  if ( sscanf(word, "%s",  label)  != 1) return;     

  word = strtok (NULL, " \t");                    /* get AXE */ 
  if ( sscanf(word, "%c", &islice) != 1) return;
  
  if (!(nslice = read_array (&vslice)) ||         /* get values and RANGE */    
      sscanf(strdup(strtok(NULL," \t]")), "%d",  &nran) != 1
      ) return;

  while (word  = strtok(NULL," \t\n"))           /* read vars list */
  {
    if (*word == '\\') {    /***** multi-line list of variables ******/
      while (fgets(line, TITLE_TIOS_LEN, tios_file) ) 
       if (*line !='%' && (word =strtok (line, " \t\n")) && *word !='%')break;
    }

    if (*word == '%')      /***** rest of the list commented out ********/
      break;
    
    if ( !(pvar = find_var_by_name(word)) ) { /*** invalid variable ******/
      printf ("TIOS-warning: not recognized variable <%s>....skipped\n", word);
      continue;
    }
    str_curr->vars.adr[str_curr->vars.count] = pvar;
    str_curr->vars.num[str_curr->vars.count++] = pvar->var_nu;
  }

  if (! str_curr->vars.count) {
    printf ("TIOS-warning: empty stream <%s>....skipped\n", label);
    return; 
  }
  str_curr->stream_nu = ++tios.str_count;
  
  pvar = str_curr->vars.adr[0]; /* make a first variable to be a stream BASE */

  pvar->strs.adr[pvar->strs.count] = str_curr;
  pvar->strs.num[pvar->strs.count++] = str_curr->stream_nu; 
  
  strcpy (str_curr->label, label);
  
  str_curr->slice.id  = islice;          /* set slice parameters */
  str_curr->slice.num = nslice;
  str_curr->slice.val = vslice;

  regrid_stream (nran, str_curr, pvar);  /* set grids for a stream or cluster*/
  set_slice (str_curr, pvar);            /* compute indexes for slices */
  
  if (str_curr->type == TP_CLUSTER) tios.clust_count++;

  if (tios.cntrl.debug) {
    printf ("Debug:%s: <%s>: %c = [ ", 
	    (str_curr->type==TP_STREAM?"Stream":"Cluster"), 
	    str_curr->label, str_curr->slice.id);
    for (i = 0; i < str_curr->slice.num; i++) 
      printf ("%.2f ", str_curr->slice.val[i]);
    printf ("], Grid=%d, %s=%d, Vars=%d\n", 
	    pvar->grid->grid_nu,(str_curr->type==TP_STREAM?"Range":"Cluster"), 
	    nran, str_curr->vars.count);
  }
  str_curr->addr = str_curr->curr = (ADDR *) calloc (1, sizeof (struct ADDR));  
  str_curr = (str_curr->next = (STREAM *) calloc (1, sizeof (struct STREAM)));
}

/*************************************************************************/
void tios_init_(tios_f, data_f, alpha, beta, gamma)
/*************************************************************************/
char *tios_f, *data_f;
REAL *alpha, *beta, *gamma;
{
  char *name, space[100];
  
  sz_REAL = sizeof(REAL);
  sz_UINT = sizeof(UINT);

  if (!(tios_file = fopen (strtok(tios_f, " \t\n"), "r")))
    {
      printf ("Cannot open file <%s>\n", tios_f);
      exit (-1);
    }

  name = strtok(data_f, " \t\n");
  strncpy (space, name, strlen(name));

  strcpy (space+strlen(name), ".data");
  if (tios.cntrl.restart) {
    if (!fopen(space, "r")) {
        printf ("Cannot open file <%s> for continuation run.\n", space);
        exit (-1); } }
  data_file = (tios.cntrl.restart ? fopen(space, "r+") : fopen(space, "w+"));

  strcpy (space+strlen(name), ".indx");
  indx_file = (tios.cntrl.restart ? fopen(space, "r+") : fopen(space, "w+"));

  tios.version = Version;
  tios.time_resolution = ENSO_DAY;
  tios.pole_alpha = *alpha ;
  tios.pole_beta =  *beta;
  tios.pole_gamma =  *gamma;
  
  tios.grids   = grid_curr = (GRID   *)calloc (1, sizeof(struct GRID)); 
  tios.maps    = map_curr  = (MAP    *)calloc (1, sizeof(struct MAP)); 
  tios.ranges  = rang_curr = (RANGE  *)calloc (1, sizeof(struct RANGE)); 
  tios.gr_maps             = (GR_MAP *)calloc (1, sizeof(struct GR_MAP)); 
  tios.vars    = var_curr  = (VAR    *)calloc (1, sizeof(struct VAR)); 
  tios.streams = str_curr  = (STREAM *)calloc (1, sizeof(struct STREAM)); 

  if (tios.cntrl.restart) {
    restart_tios(indx_file);
    fseek (indx_file, (size_t)tios.addr_start, SEEK_SET);
    fseek (data_file, (size_t)tios.addr_end, SEEK_SET);
  }
  else {
    lwrite_tios(data_file);
    lwrite_tios(indx_file);
    tios.addr_start = tios.addr_grids = ftell(indx_file);
  }
  tios.cntrl.updated = 0;
}

/*************************************************************************/
void tios_map_(addr, NXY, NMAP, comp)
/*************************************************************************/
int *addr, *NXY, *NMAP, *comp;
{
  int i;

  map_curr->map_nu = ++tios.map_count;

  map_curr->NXY  = *NXY;
  map_curr->NMAP = *NMAP;

  map_curr->cmp = (int *)calloc ((size_t)*NXY, sizeof(int)); 
/* this routine currently builds the following map (others can be easily added)
   comp(1:compressed) - FORTRAN-like indexes of valid points in XY chunk,
   comp(i) !=0 for all (i).
   map->cmp(1:full)   - fill array of FORTRAN-like indexes of compressed list
   cmp(i) == 0 for a missing data point; !=0 for an index of a compressed array
*/
  for (i = 0; i < *NMAP; i++) map_curr->cmp[comp[i]-1] = i+1;

  *addr = map_curr->map_nu;
  map_curr = (map_curr->next = (MAP *)calloc(1, sizeof(struct MAP))); 
}

/*************************************************************************/
void tios_grid_(addr, NX, NY, NZ, xx, yy, zz)
/*************************************************************************/
int *addr, *NX, *NY, *NZ;
REAL *xx, *yy, *zz;
{
  grid_curr->grid_nu = ++tios.grid_count;

  grid_curr->NX = *NX;
  grid_curr->NY = *NY;
  grid_curr->NZ = *NZ;

  grid_curr->xx = xx;
  grid_curr->yy = yy;
  grid_curr->zz = zz;

  if (tios.cntrl.debug) {    
    printf ("Debug: Grid #%d: \n", grid_curr->grid_nu);
    printf ("       NX = %3d", *NX);
    if (*NX) 
      printf ("\txx : %.4g %.4g %.4g ... %.4g\n", 
	      xx[0],xx[1],xx[2],xx[*NX-1]);
    printf ("       NY = %3d", *NY);
    if (*NY) 
      printf ("\tyy : %.4g %.4g %.4g ... %.4g\n", 
	      yy[0],yy[1],yy[2],yy[*NY-1]);
    printf ("       NZ = %3d", *NZ);
    if (*NZ) 
      printf ("\tzz : %.4g %.4g %.4g ... %.4g\n", 
	      zz[0],zz[1],zz[2],zz[*NZ-1]);
  }
  *addr = grid_curr->grid_nu;
  grid_curr = (grid_curr->next = (GRID *)calloc (1, sizeof (struct GRID))); 
}

/*************************************************************************/
void tios_var_(var, label, igrid, imap)
/*************************************************************************/
REAL *var;
char *label;
int *igrid, *imap;
{
  GRID *pgrid;
  MAP  *pmap;
  char str[100] = "j";

  var_curr->var_nu = ++tios.var_count;
  var_curr->var = var;

  strcpy(str,label);
  strcpy (var_curr->label, strtok(str," \t\n"));

  pgrid = tios.grids;
  while (pgrid->next)
    if (pgrid->grid_nu == *igrid) {
      var_curr->grid = pgrid;
      break;
    }
    else
      pgrid = pgrid->next;

  if (!pgrid) {
    printf("TIOS:tios_var: Unknown grid for variable <%s>\n", label);
    exit (-1);
  }
  
  if (tios.cntrl.debug) printf ("Debug: set variable <%s> on the grid %d\n",
			  var_curr->label, var_curr->grid->grid_nu);
  if (*imap) {
    pmap = tios.maps;
    while (pmap->next)
      if (pmap->map_nu == *imap) {
	var_curr->map = pmap;
	break;
      }
      else
	pmap = pmap->next;

    if (!pmap) {
      printf("TIOS:tios_var: Unknown map for variable <%s>\n", label);
      exit (-1);
    }
    var_curr->flag = (REAL)AUTO_FLAG;
  }

  var_curr = (var_curr->next = (VAR *) calloc (1, sizeof (struct VAR)));
}

/*************************************************************************/
int tios_idvar_(label, igrid, imap)
/*************************************************************************/
char *label;
int *igrid, *imap;
{
  GRID *pgrid;
  MAP  *pmap;
  char str[100] = "j";

  var_curr->var_nu = ++tios.var_count;
  
  strcpy(str,label);

  strcpy (var_curr->label, strtok(str," \t\n"));

  pgrid = tios.grids;
  while (pgrid->next)
    if (pgrid->grid_nu == *igrid) {
      var_curr->grid = pgrid;
      break;
    }
    else
      pgrid = pgrid->next;

  if (!pgrid) {
    printf("TIOS:tios_idvar: Unknown grid for variable <%s>\n", label);
    exit (-1);
  }

  if (tios.cntrl.debug) printf ("Debug: set variable <%s> on the grid %d\n",
			  var_curr->label, var_curr->grid->grid_nu);
  if (*imap) {
    pmap = tios.maps;
    while (pmap->next)
      if (pmap->map_nu == *imap) {
	var_curr->map = pmap;
	break;
      }
      else
	pmap = pmap->next;

    if (!pmap) {
      printf("TIOS:tios_idvar: Unknown map for variable <%s>\n", label);
      exit (-1);
    }
    var_curr->flag = (REAL)AUTO_FLAG;
  }

  var_curr = (var_curr->next = (VAR *) calloc (1, sizeof (struct VAR)));
  return (int)tios.var_count;
}

/*************************************************************************/
void tios_read_()
/*************************************************************************/
{
  VAR    *pvar;
  GRID   *pgrid;
  GR_MAP *pgm;
  STREAM *pstr;
  char line0[TITLE_TIOS_LEN], *word;
  static char *line;
  static int  lenline;

  while (fgets(line0, TITLE_TIOS_LEN, tios_file)) {

    if (lenline < strlen(line0)) lenline = strlen(line = strdup(line0)); 
    else                         strcpy(line, line0);

    if (*line == '%' || !(word = strtok (line, " \t")))
      continue;
    
    else if (!strncmp (word, "DEFINE", strlen("DEFINE"))) 
      set_define ();
    
    else if (!strncmp (word, "RANGE", strlen("RANGE"))) 
      set_range (line0, word);

    else if (!strncmp (word, "CLUSTER", strlen("CLUSTER"))) 
      set_cluster (line0);
    
    else if (!strncmp (word, "STREAM", strlen("STREAM"))) 
      set_stream (line0);
  }

  if (tios.cntrl.restart) {
    fseek (indx_file, (size_t)tios.addr_strs, SEEK_SET);
    pstr = tios.streams;
    while (pstr->next) {
      restart_stream (pstr);
      pstr = pstr->next;
    }
  }
  else {
                     /* write grids, maps, vars structures into indx_file  */
    fseek (indx_file, (size_t)tios.addr_grids, SEEK_SET);
    pgrid = tios.grids;
    while (pgrid->next) {
      lwrite_grid (pgrid);
      pgrid = pgrid->next;
    }

    pgm = tios.gr_maps;
    while (pgm->next) {
      lwrite_gmap (pgm);
      pgm = pgm->next;
    }
    tios.addr_vars = ftell(indx_file);

    pvar = tios.vars;
    while (pvar->next) {
      lwrite_var (pvar);
      pvar = pvar->next;
    }
    tios.addr_strs = ftell(indx_file);
    fflush (indx_file);
    fflush (data_file);
  }
}
/*************************************************************************/
void tios_cntrl_(key, val)
/*************************************************************************/
int *key; 
void *val;
{
  int   *pi;
  float *pf;

  switch (*key) {
  case 1: /* everystep */
    pi = (int *)val;
    tios.cntrl.everystep = *pi;
    break;
  case 2: /* debug */
    pi = (int *)val;
    tios.cntrl.debug = *pi;
    break;
  case 3: /* restart */
    pi = (int *)val;
    tios.cntrl.restart = *pi;
    break;
  case 4: /* resolution */
    pf = (float *)val;
    tios.time_resolution = *pf;
    break;
  }
}

/*.......................................................................*/
static void str_time_init (ps, time)
/*.......................................................................*/
STREAM *ps;
REAL *time;
{
  REAL time0 = (REAL)0.;

  ps->time.init = 1;
  if ( ! ps->time.absolute ) time0 = *time;

  if (ps->time.step == AUTO_FLAG) ps->time.step = 1.;

  switch (ps->time.fmt) {

  case 1:                              /* format in days */
    ps->time.step *= ENSO_DAY; 
    if (ps->time.sta == AUTO_FLAG) 
      ps->time.sta = ps->time.next = *time;
    else 
      if (tios.cntrl.restart) 
	ps->time.next += ps->time.step;
      else
	ps->time.sta = ps->time.next = time0 + ps->time.sta * ENSO_DAY;
    
    ps->time.end = ((ps->time.end == AUTO_FLAG) ?MAXFLOAT: 
		     time0 + ps->time.end * ENSO_DAY);
    break;

  case 2:                              /* format in months */
    if (ps->time.sta == AUTO_FLAG) 
      ps->time.sta = ps->time.next = *time;
    else 
      if (tios.cntrl.restart) 
	ps->time.next += ps->time.step;
      else
	ps->time.sta = ps->time.next = time0 + ps->time.sta;
    
    ps->time.end = ((ps->time.end == AUTO_FLAG) ?MAXFLOAT: time0+ps->time.end);
    break;
  case 3:                              /* format in years */
    ps->time.step *= 12.;
    if (ps->time.sta == AUTO_FLAG) 
      ps->time.sta = ps->time.next = *time;
    else 
      if (tios.cntrl.restart) 
	ps->time.next += ps->time.step;
      else
	ps->time.sta = ps->time.next = time0 + 12.* ps->time.sta; 
    
    ps->time.end = ((ps->time.end == AUTO_FLAG) ? MAXFLOAT : 
		                                  time0 + 12.* ps->time.end);
    break;
  }
}  

/*.......................................................................*/
static int isit_time (time, pstr)
/*.......................................................................*/
REAL *time;
STREAM *pstr;
{
  if ( APPROX_EQ(*time, pstr->time.next, tios.time_resolution) ) {
    if (*time > pstr->time.end) 
      return 0;
    else {
      pstr->time.next += pstr->time.step; 
      if (pstr->time.next > pstr->time.end) pstr->time.next = 0.;
      return 1;
    }
  }
  else
    return 0;
}

/*-----------------------------------------------------------------------*/
int sio_putvar(pvar, time, func, action)
/*-----------------------------------------------------------------------*/
VAR  *pvar;
REAL *time;
int  action;
void (*func)();
{
  register int i;
  STREAM *pstr;
  int iret = 0;
  
  if ( !tios.time_end && !tios.cntrl.restart) tios.time_begin = *time;
  tios.time_end = *time;

  for (i = 0; i < pvar->strs.count; i++) {
    pstr = pvar->strs.adr[i];
    
    if ( ! pstr->time.init ) str_time_init (pstr, time); 
    
    if (tios.cntrl.everystep || isit_time (time, pstr)) {
      if (action) { 
	func();
	action = 0;
      }
      write_stream_data (pstr, time);
      iret = tios.cntrl.updated = 1;
    }
  }

  return iret;
}

/*************************************************************************/
int tios_putvar_(var, time, func)
/*************************************************************************/
REAL *var, *time;
void *func;
{
  VAR *pvar = tios.vars;
  int action = *((int *)func);
  int iret = 0;

  while (pvar->var != var) 
    if (!pvar->next) 
      return iret;
    else 
      pvar = pvar->next;

  return sio_putvar(pvar, time, func, action);
}

/*************************************************************************/
int tios_putidvar_(id, var, time, func)
/*************************************************************************/
int *id;
REAL *var, *time;
void *func;
{
  VAR *pvar = tios.vars;
  int action = *((int *)func), iret = 0;

  while (pvar->var_nu != *id) 
    if (!pvar->next) 
      return iret;
    else 
      pvar = pvar->next;

  pvar->var = var;
  iret = sio_putvar(pvar, time, func, action);
  pvar->var = NULL;
  return iret;
}

/*************************************************************************/
void tios_save_()
/*************************************************************************/
{
  STREAM  *pstr;
  ADDR    *paddr;
  UINT    *pint;
  REAL    *pflt;
  register int i;

  if (tios.cntrl.updated) {

    fseek (indx_file, (size_t)tios.addr_strs, SEEK_SET);

    pstr = tios.streams;
    while (pstr->next) {
      lwrite_stream (pstr);
      
      if ( pstr->time.count) {
	pint = (UINT *)malloc(sizeof(UINT) * pstr->time.count);
	pflt = (REAL *)malloc(sizeof(REAL) * pstr->time.count);

	paddr = pstr->addr;
	for (i = 0; i < pstr->time.count; i++) {
	  pint[i] = paddr->val;
	  pflt[i] = paddr->time;
	  paddr   = paddr->next;      
	}
	
	loc_wr (sz_UINT, pstr->time.count, (void *)pint, indx_file);
	loc_wr (sz_REAL, pstr->time.count, (void *)pflt, indx_file);

	fflush (indx_file);

	free(pint);
	free(pflt);
      }
      pstr = pstr->next;
    }

    tios.addr_end = ftell(data_file);
    fseek (indx_file, 0L, SEEK_SET);
    lwrite_tios (indx_file);
    fflush (indx_file);
    fflush (data_file);
    
    tios.cntrl.updated = 0;
  }
}
/*************************************************************************/
void tios_close_()
/*************************************************************************/
{
  UINT *temp;
  VAR  *pvar;
  STREAM  *pstr;
  ADDR *paddr;
  UINT *pint;
  REAL *pflt;
  register int i;

  if (!tios.cntrl.updated) return;

  pvar = tios.vars;
  while (pvar->next) {
    temp = (UINT *)pvar;
    pvar = pvar->next;
    free ((VAR *)temp);
  }
  
  fseek (indx_file, (size_t)tios.addr_strs, SEEK_SET);
  pstr = tios.streams;
  while (pstr->next) {
    lwrite_stream (pstr);
    
    if (pstr->time.count) {
      pint = (UINT *)malloc(sz_UINT*pstr->time.count);
      pflt = (REAL *)malloc(sz_REAL*pstr->time.count);
      
      for (paddr = pstr->addr,i = 0; i < pstr->time.count; i++) {
	pint[i] = paddr->val;
	pflt[i] = paddr->time;
	temp = (UINT *)paddr;
	paddr = paddr->next;
	free ((ADDR *)temp);
      }
      
      loc_wr (sz_UINT, pstr->time.count, (void *)pint, indx_file);
      loc_wr (sz_REAL, pstr->time.count, (void *)pflt, indx_file);

      fflush (indx_file);

      free(pint);
      free(pflt);
    }

    temp = (UINT *)pstr;
    pstr = pstr->next;
    free ((STREAM *)temp);
  }

  tios.addr_end = ftell(data_file);
  fseek (indx_file, 0L, SEEK_SET);
  lwrite_tios (indx_file);

  close (data_file);
  close (indx_file);
}
/******************************************************************/
