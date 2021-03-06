# https://opensource.com/article/18/8/what-how-makefile#:~:text=The%20make%20utility%20requires%20a,be%20installed%20using%20make%20install%20.
ARCH = Linux
FC    = x86_64-linux-gnu-gfortran-4.8 #gfortran  # works fine for gfortran-4.8 or gfortran-5.5
CC    = x86_64-linux-gnu-gcc-4.8 # gcc
LCDF   = -lnetcdf -lnetcdff

ifeq (, $(shell which x86_64-linux-gnu-gfortran-4.8))
FC = gfortran
endif

ifeq (, $(shell which x86_64-linux-gnu-gcc-4.8))
CC = gcc
endif

# gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44)
# GNU Fortran (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44)

FFLAGS = -ffixed-line-length-132 -Ofast
CFLAGS = -Wno-implicit-function-declaration -Ofast

MODEL  = tcom

default: $(MODEL)

OBJ0    = om_main.o
OBJM  = wrap-mod.o data-mod.o sst-mod.o
OBJF  = om_core.o om_ekm.o om_sst.o om_leap.o \
om_equi.o om_forc.o om_qflux.o om_tios.o om_mem.o om_wrap.o fodb.o 
OBJC    = om_c.o codb.o daio.o sio.o

help:
	@echo '.  Use "make" in order to compile:\n'

LDOPT   = 

LIBS    = $(LCDF)

%.o: %.F
	$(FC) $(FFLAGS) -c $*.F

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c

$(MODEL):	$(OBJM) $(OBJF) $(OBJ0) $(OBJC)
	$(FC) $(FFLAGS) $(LDOPT) -o $(MODEL) om_main.o $(OBJM) $(OBJF) $(OBJC) $(LIBS)

tios2cdf: fodb.o $(OBJC)
	$(CC) $(CFLAGS) -o tios2cdf tios2cdf.c fodb.o $(OBJC) $(LIBS) -lgfortran

all: $(MODEL) tios2cdf

om_forc.o: wdata.mod wsst.mod wrap.mod
om_core.o om_forc.o om_sst.o om_tios.o : om_core.h 
om_forc.o om_mem.o om_wrap.o           : om_data.h
om_forc.o om_main.o om_mem.o om_tios.o om_wrap.o: om_wrap.h
om_forc.o om_mem.o om_wrap.o           : data-mod.o
om_forc.o om_main.o om_mem.o om_tios.o om_wrap.o: wrap-mod.o

om_equi.o                             : om_equi.h
om_main.o om_core.o om_wrap.o         : om_para.h
om_ekm.o om_main.o om_sst.o om_tios.o : om_sst.h
om_ekm.o om_main.o om_sst.o om_tios.o : sst-mod.o

#............................................end of Makefile
