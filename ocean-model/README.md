# Ocean model.

The ocean model code is built on legacy Fortran 90 and C code.

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

```
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