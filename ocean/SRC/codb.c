#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>

void odb_creat_ 
(name, ln, key)
char *name;
int  *ln, *key;
/******************************************************************* 
   An addition to HDF-NetCDF library:
   if key = 0, will create an empty NetCDF file if the one which is 
               specified doesn't exist.
   if key = 1 will create an empty NetCDF file unconditionally.
   ln is the length of the path "name".

   Senya Basin, 1993.
********************************************************************/
{
  char *pb = strdup(name), *pn = pb;
  int  file, count = 28;

  while (isprint(*pn) && (!isspace(*pn)) && pn < pb + *ln) pn++; 
  *pn = '\0';

  if (*key) { 
    if ((file = open(pb, O_RDWR | O_CREAT | O_TRUNC, 0644)) == -1) {
      printf ("odb_creat: can't create <%s>.\n", pb);
      exit (-1);
    }
  }
  else 
    if ((file = open(pb, O_RDWR | O_CREAT | O_EXCL, 0644)) == -1) return;

  write (file, "CDF\001", 4L);
  while (count--) write (file, "\000", 1L);

  close(file);
}













