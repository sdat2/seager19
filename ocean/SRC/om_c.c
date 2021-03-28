/**************************************************************************
*  This is a library to provide an *easy* 
*  ASCII formatted input for FORTRAN/C.
*
*  Senya Basin, 1994
**************************************************************************
$Source: /ghome/cvsroot/work/tios/dyn_c.c,v $
$Author: naomi $
$Revision: 1.1.1.1 $
$Date: 2003/08/12 16:01:47 $
$State: Exp $
***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#define MAX_LINE    200
/******* FORTRAN-word to C-string pointer *******/
#define fw_to_cp(p) strtok(strdup(p)," ,;\t\n\r\f\v\b\000") 
#define get_word(p) strtok(p," \t[](){},;\n\r\f\v\b\000") 
#define get_arra(p) strtok(p," \t,;[({\n\r\f\v\b\000") 

static FILE *file, *ftrace;
static char buff1[MAX_LINE], buff2[MAX_LINE], *psection;
static long psec_pos;
static int  trace_input, search_in_main;

static int found(int usends, long offset, char *tag)
/*******************************************************/
{ char *word;

  fseek (file, offset, SEEK_SET);

  while (fgets(buff1, MAX_LINE, file)) {
    
    if (strchr("%#",*buff1) || !(word = get_word(strcpy(buff2, buff1))) ) 
      continue;
    else if (usends && !strcasecmp(word, "end")) 
      return 0;
    else if (*word == '+' && !strcmp(word+1, tag)) 
      return 1;
  }
  return 0;
}

static int tag_found (char *tag)
/********************************************************/
{
  return ( (psec_pos && found(1, psec_pos, tag)) || 
	   (search_in_main && found(1, 0L, tag)) );
}

void prtstop1 (mess, s1)
char *mess, *s1;
{
  fprintf (stderr, mess, s1);
  exit(-1);
}

void inp_file_(name)
/************* Set the name of input FILE *********************/
char *name;
{
  if (file) fclose(file);
  if ( !(file = fopen(fw_to_cp(name), "r"))) 
    prtstop1 ("!!!inp_file: can't open <%s> for reading\n",fw_to_cp(name));
  psection = NULL;
  psec_pos = 0L;
  search_in_main = 1;
  trace_input = 0;
}

void inp_trace_(name)
/************* Set the name of trace file for output **********/
char *name;
{
  if (ftrace) fclose(ftrace);
  if ( !(ftrace = fopen(fw_to_cp(name), "w+"))) 
    prtstop1 ("!!!inp_trace: can't open <%s>.\n",fw_to_cp(name));
  trace_input = 1;
}

void inp_vrnt_(name, num)
/************* Set the VARIANT as a first place to search **/
char *name;
int  *num;
{
  char *str, space[MAX_LINE];

  str = fw_to_cp(name);
  strncpy (space, str, strlen(str));

  if (*num >= 0) sprintf (space+strlen(space), "_%1u", *num);

  psection = strdup(space);
  if (found(0,0L,psection)) psec_pos = ftell(file);
  else                      psec_pos = 0L;
  search_in_main = 1;
}

void inp_sect_(name)
/************* Set the SECTION as a first place to search **/
char *name;
{
  char *str, space[MAX_LINE];

  str = fw_to_cp(name);
  strncpy (space, str, strlen(str));

  psection = strdup(space);
  if (found(0,0L,psection)) psec_pos = ftell(file);
  else                      psec_pos = 0L;
  search_in_main = 0;
}

int inp_int_(tag, dflt)
/*********** input an INTEGER number ***********************/
char *tag;
int *dflt;
{ int val;

  if (tag_found(fw_to_cp(tag))) {
    if ( sscanf(get_word(NULL), "%d", &val) != 1) 
      prtstop1 ("!!!inp_int: can't read <%s>\n", fw_to_cp(tag));
  }
  else
    val = *dflt;

  if (trace_input) 
    fprintf(ftrace, "+%-20s %d\n", fw_to_cp(tag), val), fflush(ftrace);
  return val;
}

float inp_flt_(tag, dflt)
/************* input a REAL/FLOAT number *******************/
char *tag;
float *dflt;
{ float val;

  if (tag_found(fw_to_cp(tag))) {
    if ( sscanf (get_word(NULL), "%g", &val) != 1) 
      prtstop1 ("!!!inp_flt: can't read <%s>\n", fw_to_cp(tag));
  }
  else 
    val = *dflt;

  if (trace_input) 
    fprintf(ftrace, "+%-20s %g\n", fw_to_cp(tag), val), fflush(ftrace);
  return val;
}

double inp_dbl_(tag, dflt)
/************** input a REAL*8/DOUBLE number ****************/
char   *tag;
double *dflt;
{ double val;

  if (tag_found(fw_to_cp(tag))) {
    if ( sscanf (get_word(NULL), "%lg", &val) != 1) 
      prtstop1 ("!!!inp_flt: can't read <%s>\n", fw_to_cp(tag));
  }
  else 
    val = *dflt;

  if (trace_input) 
    fprintf(ftrace, "+%-20s %lg\n", fw_to_cp(tag), val), fflush(ftrace);
  return val;
}

int inp_str_(tag, dflt, val)
/************  input a 'CHARACTER' "string" *****************/
char *tag, *dflt, *val;
{ char *p1, *p2;

  if (tag_found(fw_to_cp(tag))) {
    
    if (!(p1 = strpbrk (buff1, "\"\'")) || !(p2 = strchr(p1+1,*p1)) )
      prtstop1 ("!!!inp_str: error reading <%s>\n", fw_to_cp(tag));
    
    *p2 = '\0';
    strcpy(val, &p1[1]);
  }
  else 
    if (val != dflt) strcpy(val, dflt);

  if (trace_input) 
    fprintf(ftrace, "+%-20s \"%s\"\n", fw_to_cp(tag), val), fflush(ftrace);
  return strlen(val);
}

int inp_iarr_ (tag, ddim, darr, arr)
/************* input an integer ARRAY ************************/
char *tag;
int  *ddim, *darr, *arr;
{ register int i; char *word; int dim;
  
  if (tag_found(fw_to_cp(tag))) {
    dim = 0;
    while ((word = get_arra(NULL)) || 
	   (fgets(buff1,MAX_LINE,file) && !strchr("%#",*buff1) &&
	    (word = get_arra(strcpy(buff2,buff1))) )
	  ) 
      {
	if ( strpbrk(word, ")]}")) {
	  sscanf(word, "%d", &arr[dim++]);
	  break;
	}
	if ( sscanf(word, "%d", &arr[dim++]) != 1) 
	  prtstop1 ("!!!inp_iarr: can't read <%s>\n", fw_to_cp(tag));
      }
  }
  else {
    dim = *ddim;
    if (arr != darr) for(i=0; i<dim; i++) arr[i]=darr[i];  
  }

  if (trace_input) {
    fprintf(ftrace, "+%-20s [", fw_to_cp(tag));
    for(i = 0; i<dim;) {
      fprintf(ftrace, "%d ", arr[i]);
      if (!(++i%10)) fprintf(ftrace,"\n%23c",' ');
    }
    fprintf(ftrace, "]\n"), fflush(ftrace);
  }
  return dim;
}

int inp_rarr_ (tag, ddim, darr, arr)
/************* input an REAL ARRAY ************************/
char  *tag;
int   *ddim;
float *darr, *arr;
{ register int i; char *word,*pntr; int dim;
  
  if (tag_found(fw_to_cp(tag))) {
    dim = 0;
    while ((word = get_arra(NULL)) || 
	   (fgets(buff1,MAX_LINE,file) && !strchr("%#",*buff1) &&
	    (word = get_arra(strcpy(buff2,buff1))) )
	  ) 
      {
	if ( pntr = strpbrk(word, ")]}") ) {
	  if (pntr != word) sscanf(word, "%g", &arr[dim++]);
	  break;
	}
	if ( sscanf(word, "%g", &arr[dim++]) != 1) 
	  prtstop1 ("!!!inp_rarr: can't read <%s>\n", fw_to_cp(tag));
      }
  }
  else {
    dim = *ddim;
    if (arr != darr) for(i=0; i<dim; i++) arr[i]=darr[i];  
  }

  if (trace_input) {
    fprintf(ftrace, "+%-20s [", fw_to_cp(tag));
    for(i = 0; i < dim;) {
      fprintf(ftrace, "%g ", arr[i]);
      if (!(++i%10)) fprintf(ftrace,"\n%23c",' ');
    }
    fprintf(ftrace, "]\n"), fflush(ftrace);
  }
  return dim;
}

int inp_sarr_ (tag, ddim, darr, dlen, alen, arr)
/************* input a STRING ARRAY ************************/
char  *tag;
int   *dlen, *ddim, *alen;
char  *darr, *arr;
{ register int i; int len, dim, slen;
  char *p1, *p2;

  if (tag_found(fw_to_cp(tag))) {
    dim = 0; p1 = buff1;
    while (((p1 = strpbrk(p1, "\"\'")) && (p2 = strchr(p1+1,*p1)) ) || 
	   ((p1 = fgets(buff1,MAX_LINE,file)) && !strchr("%#", *p1) &&
           ((p1 = strpbrk(p1, "\"\'")) && (p2 = strchr(p1+1,*p1)))) )
      {
        for (i=0; i<*dlen; i++) arr[*dlen*dim+i] = '\0';
	slen = (int)(p2-p1-1);
	strncpy(&arr[*dlen*dim], &p1[1], (*dlen>slen ? slen : *dlen));
	alen[dim] = strlen(&arr[*dlen*dim]);
        dim++; p1 = p2+1;
      }
  }
  else {
    dim = *ddim;
    if (arr != darr) for(i=0; i<*dlen*dim; i += *dlen)  
      strncpy(&arr[i],&darr[i],*dlen);
    for(i = 0; i < dim; i++)  {
/*      alen[i] = ((len = strlen(&arr[*dlen*i])) > *dlen) ? *dlen : len; */
      p1 = strchr(&arr[*dlen*i],' ');
      *p1 = '\0';
      alen[i] = strlen(&arr[*dlen*i]);
    }
  }

  if (trace_input) {
    fprintf(ftrace, "+%-20s [", fw_to_cp(tag));
    for(len = 23,i = 0; i < dim;) {
      slen = strlen(&arr[*dlen*i]);
      if (slen > *dlen) slen = *dlen;
      fprintf(ftrace, "\"%*s\" ", slen, &arr[*dlen*i]);
      len += 2 + slen;
      if (!(++i%10) || len > 70) len = 23,fprintf(ftrace,"\n%23c",' ');
    }
    fprintf(ftrace, "]\n"), fflush(ftrace);
  }
  return dim;
}

void inp_any_(tag, dflt, val, type)
/************ input *any* OBJECT ***********************/
char *tag, *type;
void *dflt, *val;
{
  int sz;
  char *p1, *p2, fmi[12], fmo[16], 
       *word  = strdup(fw_to_cp(type));

  strcpy (fmi, "%d");
  strcpy (fmo, "+%-20s %d");

  if (!strcasecmp(word,"i")       || 
      !strcasecmp(word,"int")     || 
      !strcasecmp(word,"integer")) sz = sizeof(int); 
  else if
     (!strcasecmp(word,"f")       || 
      !strcasecmp(word,"flt")     || 
      !strcasecmp(word,"float")   || 
      !strcasecmp(word,"r")       || 
      !strcasecmp(word,"real")) sz = sizeof(float),fmi[1] = 'g', fmo[8] = 'g';
  else if
     (!strcasecmp(word,"d")       || 
      !strcasecmp(word,"dbl")     || 
      !strcasecmp(word,"dble")    || 
      !strcasecmp(word,"double"))  
      sz = sizeof(double),strcpy(fmi,"%lg"), strcpy(fmo, "+%-20s %lg");
  else if
     (!strcasecmp(word,"c")       || 
      !strcasecmp(word,"c1")      || 
      !strcasecmp(word,"char")    || 
      !strcasecmp(word,"char1")   || 
      !strcasecmp(word,"char*1")  || 
      !strcasecmp(word,"character")) sz = 1,fmi[1] = 'c', fmo[8] = 'c';
  else if
     (!strcasecmp(word,"w")       || 
      !strcasecmp(word,"word")) sz = 0,fmi[1] = 's', fmo[8] = 's';
  else if
     (!strcasecmp(word,"s")       || 
      !strcasecmp(word,"str")     || 
      !strcasecmp(word,"string")) 
      sz = 0,fmi[1] = 's', strcpy(fmo, "+%-20s \"%s\"");
  else if
     (!strcasecmp(word,"i1")      || 
      !strcasecmp(word,"int1")    || 
      !strcasecmp(word,"integer*1")) sz = 1;
  else if
     (!strcasecmp(word,"i2")      || 
      !strcasecmp(word,"int2")    || 
      !strcasecmp(word,"integer*2")) sz = 2;
  else if
     (!strcasecmp(word,"i4")      || 
      !strcasecmp(word,"int4")    || 
      !strcasecmp(word,"integer4")|| 
      !strcasecmp(word,"integer*4")) sz = 4;
  else if
     (!strcasecmp(word,"r4")      || 
      !strcasecmp(word,"real4")   || 
      !strcasecmp(word,"real*4"))    sz = 4,fmi[1] = 'g', fmo[8] = 'g';
  else if
     (!strcasecmp(word,"r8")      || 
      !strcasecmp(word,"real8")   || 
      !strcasecmp(word,"real*8"))     
      sz = 8,strcpy(fmi,"%lg"), strcpy(fmo, "+%-20s %lg");
  else if
     (!strcasecmp(word,"h")       || 
      !strcasecmp(word,"x")       || 
      !strcasecmp(word,"hex")     || 
      !strcasecmp(word,"hexadecimal"))
      sz = 1, strcpy(fmi,"%*2c%x"), strcpy(fmo, "+%-20s %#x");
  else if
     (!strcasecmp(word,"o")       || 
      !strcasecmp(word,"oct")     || 
      !strcasecmp(word,"octal"))      
      sz = 1, strcpy(fmi,"%*c%o"), strcpy(fmo, "+%-20s %#o");
  else
    prtstop1("!!!inp_any: unknown format <%s>\n", word);

  if ( tag_found(fw_to_cp(tag)) ) {
    
    if (sz) { 
      if (sscanf (get_word(NULL), fmi, val) != 1)
	prtstop1 ("!!!inp_any: can't read <%s>\n", fw_to_cp(tag));
    }
    else {
      if (*word == 'w') {
	if (word = get_word(NULL)) 
	  strcpy(val, word);
	else 
	  prtstop1 ("!!!inp_any: can't read <%s>\n", fw_to_cp(tag));
      }
      else {
	if (!(p1 = strpbrk (buff1, "\"\'")) || !(p2 = strchr(p1+1,*p1)))
	  prtstop1 ("!!!inp_any: error reading <%s>\n", fw_to_cp(tag));
	
	*p2 = '\0';
	strcpy(val, &p1[1]);
      }
    }
  }
  else
    if (val != dflt) {
      if (sz) 
	memmove(val, dflt, (size_t)sz);
      else {
	if (*word == 'w') strcpy (val, fw_to_cp(dflt));
	else              strcpy (val, dflt);
      }
    }
  if (trace_input) 
    fprintf(ftrace, fmo, fw_to_cp(tag), val), fflush(ftrace);
}

float inp_days_(tag, dflt)
/************* input a TIME in days *******************/
char *tag;
float *dflt;
{ float old, val; char *word, *fmt = "day";

  if (tag_found(fw_to_cp(tag))) {
    if ( sscanf(get_word(NULL), "%g", &old) == 1) 
      word = get_word(NULL);
    else
      prtstop1 ("!!!inp_days: can't read <%s>\n", fw_to_cp(tag));
  }
  else
    val = old = *dflt, word = fmt;

  switch (tolower(*word) ) {
  case 'h':
    val = old/24.;
    fmt = "hour";
    break;
  case 'w':
    val = 7.*old;
    fmt = "week";
    break;
  case 'm':
    val = (365./12.)*old;
    fmt = "month";
    break;
  case 'y':
    val = 365.*old;
    fmt = "year";
    break;
  default:
    val = old;
    break;
  }

  if (trace_input)
    fprintf(ftrace, "+%-20s %g  %s%c\n", fw_to_cp(tag), old, fmt,
	    ((int)old == 1 ?' ':'s')), 
    fflush(ftrace);

  return val;
}

void inp_date_(tag, dm, dd, dy, im, id, iy)
char *tag;
int *dm, *dd, *dy, *im, *id, *iy;
{
  if (tag_found(fw_to_cp(tag))) {
    if ( sscanf(get_word(NULL), "%d", im) != 1 ||
	 sscanf(get_word(NULL), "%d", id) != 1 ||
	 sscanf(get_word(NULL), "%d", iy) != 1 ) 
      prtstop1 ("!!!inp_date: can't read <%s>\n", fw_to_cp(tag));
  }
  else
    *im = *dm, *id = *dd, *iy = *dy;

  if (trace_input) 
   fprintf(ftrace,"+%-20s %d %d %d\n",fw_to_cp(tag),*im,*id,*iy),fflush(ftrace);
}

int inp_def_ (tag) 
char *tag;
{ int val;
  val = tag_found(fw_to_cp(tag)); 
  if (trace_input) 
    fprintf(ftrace, "+%-20s %d\n", fw_to_cp(tag), val), fflush(ftrace);
  return val;
}

int inp_inxt_(dflt)
/*********** input a next INTEGER from the previous TAG ****/
int *dflt;
{ int val;

  if ( sscanf(get_word(NULL), "%d", &val) != 1) 
    prtstop1 ("!!!inp_inxt: can't read <%s>\n", buff1);
  else
    val = *dflt;
  return val;
}

float inp_fnxt_(dflt)
/*********** input a next FLOAT from the previous TAG ****/
float *dflt;
{ float val;

  if ( sscanf(get_word(NULL), "%g", &val) != 1) 
      prtstop1 ("!!!inp_fnxt: can't read <%s>\n", buff1);
  else
    val = *dflt;
  return val;
}

int inp_wnxt_ (val)
/************* read a next *word* form the file ************************/
char *val;
{ char *word;

  return ((word = get_word(NULL)) || 
	 (fgets(buff1,MAX_LINE,file) && !strchr("%#",*buff1) &&
	  (word = get_word(strcpy(buff2,buff1))) )
	 ) ? strlen(strcpy(val,word)) : 0;
}

