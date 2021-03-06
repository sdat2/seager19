# Original python atmosphere model from Naomi Henderson

There are two original notebooks and python scripts,
that which are adaptations of one another. These
were both refacotred into `src/models/atmos` so 
that they could share a common input struct and 
be called systematically as an object from 
within a coupled model run.


File structure:

```
  TCAM.ipynb
  TCAM.py
  windsFromSST-K10-eps0.75.eps
  dQ.py
  README.md
  file_struct.txt
  dQ.ipynb
  DATA/
    ps-ECMWF-clim.nc
    ts-ECMWF-clim60.nc
    sfcWind-ECMWF-clim60.nc
    pr-ECMWF-trend.nc
    sst-ECMWF-clim.nc
    mask-360x181.nc
    ts-ECMWF-clim.nc
    mask-360x180.nc
    clt-ECMWF-clim60.nc
    ts-ECMWF-trend.nc
    sst-ECMWF-trend.nc
    rh-fixed-clim60.nc
    rh-ECMWF-clim60.nc
    sfcWind-ECMWF-clim.nc
    pr-ECMWF-clim.nc
  tmp/
    S91-Hq1800-PrcpLand1.nc
    S91-Hq1800-PrcpLand0.nc
```
