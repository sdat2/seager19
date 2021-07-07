
# Initial Code/Data structure

- This section shows the initial structure of the code and
   data.

- From inspection, files with the same names in different
  dircetories are identical.

- How to couple the ocean and atmospheric model is not
   immediately apparent.

## Summary

### Code Makeup

The Python code for the atmosphere model is in a Juypter Notebook.
The ocean model code is built on legacy Fortran 90 and C code.

If cloc does not exist

```bash
  sudo apt-get install cloc
```

From running the command

```bash
  cloc $(git ls-files)
```

The initial state off the code was:

```txt
    48 text files.
    45 unique files.s
    14 files ignored.
```

github.com/AlDanial/cloc v 1.84  T=0.10 s (349.4 files/s, 149217.5 lines/s)

 | Language                |       files       |     blank      |    comment      |       code |
 | ----------------------- | ----------------- | -------------- | --------------- | ---------- |
 | Fortran 77              |          15       |      1364      |       1365      |       6170 |
 | C                       |           5       |       493      |        200      |       2746 |
 | Jupyter Notebook        |           2       |         0      |        517      |        474 |
 | Python                  |           2       |       172      |        100      |        397 |
 | C/C++ Header            |           8       |        88      |         18      |        365 |
 | make                    |           1       |        15      |          1      |         36 |
 | Markdown                |           1       |         0      |          0      |          1 |
 | SUM:                    |          34       |      2132      |       2201      |      10189 |

### Code structure

The code is structured into folders:

```txt
   |-animations
   |-atmos
   |---README.md --> lists file structure of this model.
   |---DATA
   |---tmp
   |-ocean
   |---README.md --> lists file structure of this model.
   |---DATA
   |---RUN
   |-----run-model
   |-----DATA
   |-----output
   |---SRC
   |-----DATA
   |-----output
   |-requirements
```

## Detailed

- File by file structure.

### Data

- The data is currently not stored in the github repository, as it takes up roughly 3.5 GB.
- Given the duplication of data, it should be possible to
   reduce the ammount to a managable ammount, and either
   keep it on `git-lfs` or create
   an automatic import script from Dropbox.

```txt
atmos
  DATA/      # everything ending with clim60 read by dQ.py
    ps-ECMWF-clim.nc          # read in by TCAM.py
    ts-ECMWF-clim60.nc        # read in by dQ.py
    sfcWind-ECMWF-clim60.nc   # read in by dQ.py
    pr-ECMWF-trend.nc         # read in by TCAM.py 
    sst-ECMWF-clim.nc         # read in by TCAM.py
    mask-360x181.nc
    ts-ECMWF-clim.nc          # read in by TCAM.py
    mask-360x180.nc           # mask
    clt-ECMWF-clim60.nc       # read in by dQ.py
    ts-ECMWF-trend.nc         # read in by TCAM.py
    sst-ECMWF-trend.nc        # read in by TCAM.py
    rh-fixed-clim60.nc        # read in by dQ.py
    rh-ECMWF-clim60.nc        # read in by dQ.py
    sfcWind-ECMWF-clim.nc     # read in by TCAM.py
    pr-ECMWF-clim.nc          # read in by TCAM.py
  tmp/
    S91-Hq1800-PrcpLand1.nc
    S91-Hq1800-PrcpLand0.nc
ocean
  RUN
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
    output
      om_spin.20y.restart
      om_diag.nc
      om_test.indx
      om_diag.save
      om_diag.data
      om_test.data
      om_test.save
      om_diag.indx
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
```

### Code

- Duplication of code between `jupyter-notebook`s and the `python` scripts.
- Executables `om_run2f`, `tios2cdf`.
- The mixture of `Fortran77`/`C` files in the ocean model is currently pretty inpenetrable.

```txt
atmos
  TCAM.ipynb
  TCAM.py
  dQ.py
  dQ.ipynb
  DATA/
  tmp/
ocean
  RUN
    run-model
    diag.tios
    om_diag.log
    om_run2f
    qflx.ing
    output/
    month.tios
    spin.tios
    om_diag
    om_spin
    DATA/
  DATA/
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
    output/
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
    DATA/
    tios.h
    om_diag.tr
    daio.c
    om_test.log
```
