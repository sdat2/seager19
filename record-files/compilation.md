## Current compilation:

```
gfortran  -ffixed-line-length-132 -Ofast -c om_core.F
gfortran  -ffixed-line-length-132 -Ofast -c om_ekm.F
gfortran  -ffixed-line-length-132 -Ofast -c om_sst.F
gfortran  -ffixed-line-length-132 -Ofast -c om_leap.F
gfortran  -ffixed-line-length-132 -Ofast -c om_equi.F
gfortran  -ffixed-line-length-132 -Ofast -c om_forc.F
gfortran  -ffixed-line-length-132 -Ofast -c om_qflux.F
gfortran  -ffixed-line-length-132 -Ofast -c om_tios.F
gfortran  -ffixed-line-length-132 -Ofast -c om_mem.F
gfortran  -ffixed-line-length-132 -Ofast -c om_wrap.F
gfortran  -ffixed-line-length-132 -Ofast -c fodb.F
gfortran  -ffixed-line-length-132 -Ofast -c om_main.F
gcc -Wno-implicit-function-declaration -Ofast -c om_c.c
gcc -Wno-implicit-function-declaration -Ofast -c codb.c
gcc -Wno-implicit-function-declaration -Ofast -c daio.c
gcc -Wno-implicit-function-declaration -Ofast -c sio.c
gfortran  -ffixed-line-length-132 -Ofast  -o tcom om_main.o wrap-mod.o data-mod.o sst-mod.o om_core.o om_ekm.o om_sst.o om_leap.o om_equi.o om_forc.o om_qflux.o om_tios.o om_mem.o om_wrap.o fodb.o  om_c.o codb.o daio.o sio.o -lnetcdf -lnetcdff
```

