# Ocean model

The ocean model code is built on legacy Fortran 90 and C code.

TDEEP / tdeep - temperature of deep water
HTHERM / htherm - height of the thermocline?
W1 / w1 - upwelling speed

To backup the data:

```bash
sh ./backup-data.sh
```

```bash
sh ./comp-data.sh >> changes-mid-run.txt
```

```bash
sh ./comp-dat.sh
```

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

```yaml
before_install:
  - sudo apt-get install gfortran

matrix:
  include:
    # works on Precise and Trusty
    - os: linux
      addons:
        apt:
          sources:
            - ubuntu-toolchain-r-test
            - llvm-toolchain-precise-3.6
          packages:
            - GCC-4.8.5
      env:
        - MATRIX_EVAL="CC=gcc-4.8.5 && CXX=g++-4.8.5"
```

```yaml
language: c
sudo: required
before_install:
  - sudo apt-get install gfortran

script:
  - gfortran -fprofile-arcs -ftest-coverage -O0 hello.f90 -o hello
  - ./hello

after_success:
  - bash <(curl -s https://codecov.io/bash)
```

```bash
  sudo apt-get install python-software-properties
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test
  sudo apt-get update
  sudo apt-get install gcc-4.8
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 50
  sudo apt-get install gfortran-4.8
```
