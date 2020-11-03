/*-----------------------------------------------------------------------
 This is a set of "C" /fopen/fwrite/fread/fseek/ftell/ -based 
   subroutines to provide efficient and easy direct access i/o from FORTRAN.

  ===========================================================================
   opda (lun, cntrl, file)
   clda (lun)
      - to open/close file as DIRECT ACCESS one.

   rdda (lun, len, addr, dim, array)
   wrda (lun, len, addr, dim, array)      
      - to read/write "dim" elements of the array from specified address "addr"

   integer function irdda (len, dim, array)
   integer function iwrda (len, dim, array)
      - to read/write "dim" elements of the array from current address & file.
        they both return the value of next available address.

   integer function lenda (lun) - determine size of a file in bytes. 
   integer function mfloat_read (lun, array) 
      - read ASCII formatted FLOATING POINT data, ignoring white-space 
        and punctuation characters  
   integer function mint_read (lun, array) - the same for INTEGERS
  ============================================================================
    
   integer lun        - logical unit number
   integer cntrl      - if =0: reading only; !=0: update or create 
   character*(*) file - filename
   integer len        - length of the array element in bytes
   integer addr       - on input:  the array's adress in file
                      - on output: the next available address in file
   integer dim        - dimension of the array to be read/written

   S.Basin, May 1991.
-----------------------------------------------------------------------------*/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_FILES   100
#define MAX_NAMELEN 512

static FILE *node[MAX_FILES], *curr_file;

void opda_(lun, cntrl, name)
int *cntrl, *lun;
char *name;
{
  int file;
  char cbuf[MAX_NAMELEN];
  char *cb = cbuf, *ch = name;

  while (isprint(*ch) && (!isspace(*ch)) && ch<name+MAX_NAMELEN) *cb++ = *ch++; 
  *cb = '\0';

  if (!(node[*lun]))
    if (*cntrl == 0) {
      if (!(node[*lun] = fopen(cbuf, "r"))) {
	fprintf(stderr,"opda: problems with <%s> for reading\n",cbuf);
	exit(-1);
      }
    }
    else if (*cntrl == 1) {
      if (!(node[*lun] = fopen(cbuf,"r+")) && 
	  !(node[*lun] = fopen(cbuf,"w+"))) {
	fprintf (stderr,"opda: can't open <%s>\n",cbuf);
	exit (-1);
      }
    }
    else if (*cntrl == 2) {
      if ((file = open(cbuf, O_RDWR | O_CREAT | O_TRUNC, 0644)) == -1) {
	fprintf (stderr,"opda: can't create <%s>\n",cbuf);
	exit (-1);
      }
      else
	close(file);

      if (!(node[*lun] = fopen(cbuf,"w+")) ) {
	fprintf (stderr,"opda: can't open <%s>\n",cbuf);
	exit (-1);
      }
    }
  else 
    if (!(node[*lun] = fopen(cbuf,"r+")) && 
	!(node[*lun] = fopen(cbuf,"w+"))) {
      fprintf (stderr,"opda: can't open <%s>\n",cbuf);
      exit (-1);
    }

  curr_file = node[*lun];
}

void opnda_(lun, cntrl, name, ln)
int *cntrl, *lun, *ln;
char *name;
{
  char cbuf[MAX_NAMELEN];
  if ( *ln ) {
    strncpy(cbuf, name, *ln);
    cbuf[*ln]= '\0';
    opda_(lun, cntrl, cbuf);
  }
  else {
    fprintf(stderr,"opnda: invalid filename length!\n");
    exit(-1);
  }
}

void clda_(lun)
int *lun;
{
  if (fclose (node[*lun]))
    {
      fprintf (stderr, "clda: error!\n");
      exit (-1);
    }
  node[*lun] = 0;
}

void flda_(lun)
int *lun;
{
  if (fflush (node[*lun]))
    {
      fprintf (stderr, "flda: error!\n");
      exit (-1);
    }
}

void wrda_(lun, type, addr, num, data)
int *lun, *type, *addr, *num;
void *data;
{
  fseek (node[*lun], *addr, 0L);
  if (fwrite (data, 1L, *type**num, node[*lun]) != *type**num)
    {
      fprintf (stderr, "wrda: addr=%ld error\n", *addr);
      exit (-1);
    }
  *addr = ftell(curr_file = node[*lun]);
}

void rdda_(lun, type, addr, num, data)
int *lun, *type, *addr, *num;
void *data;
{
  fseek (node[*lun], *addr, 0L);
  if (fread (data, 1L, *type**num, node[*lun]) != *type**num)
    {
      *addr = 0;
      return;
    }
  *addr = ftell(curr_file = node[*lun]);
}

int iwrda_(type, num, data)
int *type, *num;
void *data;
{
  if (fwrite (data, 1L, *type**num, curr_file) != *type**num)
    {
      fprintf (stderr, "iwrda: error\n");
      exit (-1);
    }
  return ftell(curr_file);
}

int irdda_(type, num, data)
int *type, *num;
void *data;
{
  if (fread (data, 1L, *type**num, curr_file) != *type**num)
    {
      return 0;
    }
  return ftell(curr_file);
}

int lenda_(lun)
int *lun;
{
  if (fseek (node[*lun], 0L, SEEK_END)) {
    fprintf (stderr, "lenda: file #%2d is not opened\n", *lun);
    return 0;
  }
  else
    return (int)ftell(node[*lun]); 
}

int mfloat_read_(lun, data)
int *lun;
float *data;
{
  int rc, count = 0;
  while ((rc = fscanf (node[*lun], "%g", &data[count])) == 1) count++;

  return count;
}

int mint_read_(lun, data)
int *lun;
long *data;
{
  int rc, count = 0;
  while ((rc = fscanf (node[*lun], "%ld", &data[count])) == 1) count++;

  return count;
}

int nfloat_read_(lun, nc, data)
int *lun, *nc;
float *data;
{
  int rc, cnt = 0;
  while (cnt < *nc || (rc = fscanf (node[*lun], "%g", &data[cnt]) == 1)) cnt++;

  return cnt;
}

int nint_read_(lun, nc, data)
int *lun, *nc;
long *data;
{
  int rc, cnt = 0;
  while (cnt < *nc || (rc = fscanf (node[*lun], "%ld", &data[cnt]) == 1)) cnt++;

  return cnt;
}

int killsp_(size, str)
int *size;
char *str;
{
  char *ps = str;
  while (isalnum(*ps) && (ps < str + *size)) ps++;
  if (ps < str + *size) 
    {
      *ps = '\0';
      return (int)(ps - str);
    }
  else
    return *size;
}
/*-----------------------------------------------------------------------*/



