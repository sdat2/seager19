# Ocean model

The ocean model code is built on legacy Fortran 77 and C code. The description from om_main.F is:

```fortran77
c     multimode linear equatorial ocean model, using the INC scheme
c
c        authors:     N. Henderson(Naik), with contributions from: 
c                                   B. Blumenthal (multimode)
c                                   R. Seager     (AML - advective mixed layer)
c
c        references:
c                    M. Israeli, Naik, N. and Cane, M.A., 2000. 
c             "An Unconditionally Stable Scheme for the Shallow Water Equations"
c
c                    M.B. Blumenthal and Cane, M., 1989. "Accounting for 
c              parameter uncertainties in model verification: an illustration
c              with tropical sea surface temperature," 
c              J. Phys. Oceanogr.19, 815-830.
c
c                    R. Seager, Blumenthal, M.B. and Kushnir, Y., 1995.
c             "An advective atmospheric mixed layer model for ocean
c             modeling purposes: Global simulation of surface heat fluxes",
c              J. Climate, 8, 1951-1964.
c
c                    R. Seager, Kushnir, Y. and Cane, M.A., 1995.
c             "On heat flux boundary conditions for ocean models"
c              J. Phys. Oceanogr., 25, 3219-3230.
```

## Glossary of output fields

- TDEEP / tdeep - temperature of deep water
- HTHERM / htherm - height of the thermocline?
- W1 / w1 - upwelling speed in surface layer

SST W1 UDTDX VDTDY UP_FLUX QPRIME QFC QNET

## Helpful background information

Baroclinic Rossby waves:

<https://www.youtube.com/watch?v=ycs2AbC44EU>

Barotropic Rossby waves:

<https://youtu.be/pwV54L-NXzM>

Equatorial waves:

<https://youtu.be/Tdi7lulinRg>

## Comparison and analysis of model

### To backup the data

```bash
sh ./backup-data.sh
```

### To compare the data between different folders

```bash
sh ./comp-data.sh >> changes-mid-run.txt
```

## Structure of Ocean model

Compile the programs in `SRC`.

Run the programs in `RUN`.

Ocean model file structure:

```txt
  files.txt
  README.md
  file_struct.txt
  RUN
    run-model
    diag.tios
    om_diag.log
    om_run2f
    qflx.ing
    output
      om_spin.20y.restart
      om_run2f.nc
      om_spin.nc
      om_spin.save
      om_run2f.save
      om_diag.2y.restart
      om_diag.nc
      om_diag.save
      om_diag.data
      om_diag.indx
    month.tios
    spin.tios
    om_diag
    om_spin
    DATA
      rzk.pro
      spline_ECMWF.txt
      dQdf-sample.nc
      om_mask.nc
      qflx.nc
      tau-ECMWF.y
      sst-ECMWF-clim.nc
      tau-ECMWF.x
      dQdT-sample.nc
      qflx-0.nc
      tau-ECMWF-clim.x
      tau-ECMWF-clim.y
    om_diag.tr
  DATA
    rzk.pro
    spline_ECMWF.txt
    dQdf-sample.nc
    om_mask.nc
    qflx.nc
    tau-ECMWF.y
    sst-ECMWF-clim.nc
    tau-ECMWF.x
    dQdT-sample.nc
    qflx-0.nc
    tau-ECMWF-clim.x
    tau-ECMWF-clim.y
  SRC
    om_mem.F
    wrap-mod.F
    om_data.h
    diag.tios
    om_test.tr
    om_diag.log
    om_qflux.F
    om_forc.F
    om_leap.F
    Makefile
    sst-mod.F
    netcdf.inc
    qflx.ing
    output
      om_spin.20y.restart
      om_diag.nc
      om_test.indx
      om_diag.save
      om_diag.data
      om_test.data
      om_test.save
      om_diag.indx
    cuf.h
    om_main.F
    om_core.F
    om_core.h
    om_sst.h
    README
    om_sst.F
    sio.c
    om_equi.h
    om_equi.F
    codb.c
    om_tios.F
    om_wrap.h
    tios2cdf.c
    om_wrap.F
    om_test
    om_ekm.F
    om_c.c
    fodb.F
    om_diag
    data-mod.F
    om_para.h
    DATA
      rzk.pro
      spline_ECMWF.txt
      dQdf-sample.nc
      om_mask.nc
      qflx.nc
      tau-ECMWF.y
      sst-ECMWF-clim.nc
      tau-ECMWF.x
      dQdT-sample.nc
      qflx-0.nc
      tau-ECMWF-clim.x
      tau-ECMWF-clim.y
    tios.h
    om_diag.tr
    daio.c
    om_test.log
```

### Full current compiler details

```txt
GNU Fortran (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44)
Copyright (C) 2015 Free Software Foundation, Inc.

GNU Fortran comes with NO WARRANTY, to the extent permitted by law.
You may redistribute copies of GNU Fortran
under the terms of the GNU General Public License.
For more information about these matters, see the file named COPYING
```

```txt
gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44)
Copyright (C) 2015 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```
