## Model from the Methods Appendix

### Atmos and Ocean reanalysis

- ECMWF ERA-40 1958-78 & ERA-Interim 1979-2017
  - Wind at 2m from surface, 
  - precipitation, 
- ORAS4 1958-2017
  - SST, 
  - themocline depth (20C isotherm).

![gifs/clt_in_clt-ECMWF-clim60.gif](gifs/clt_in_clt-ECMWF-clim60.gif)

![gifs/dQdf_in_dQdf-sample.gif](gifs/dQdf_in_dQdf-sample.gif)

![gifs/pr_in_pr-ECMWF-clim.gif](gifs/pr_in_pr-ECMWF-clim.gif)

![gifs/ps_in_ps-ECMWF-clim.gif](gifs/ps_in_ps-ECMWF-clim.gif)

![gifs/rh_in_rh-ECMWF-clim60.gif](gifs/rh_in_rh-ECMWF-clim60.gif)

![gifs/sst_in_sst-ECMWF-clim.gif](gifs/sst_in_sst-ECMWF-clim.gif)

![gifs/ts_in_ts-ECMWF-clim60.gif](gifs/ts_in_ts-ECMWF-clim60.gif)


### Atmospheric Model

  (u, v, w) = (u'cos(pi z/zT), v'cos(pi z/zT), w'sin(pi z/zT))

z is the height, and zT is the top of the troposphere.

  (theta, Q) = (theta', Q')(theta0/theta00) sin(pi z/zT)

  p = p'(rho0 / rho00) cos(pi z/zT)

  (p/rho0)_z = g theta / theta00

  p'/rho00 = (g zT/ pi theta00)theta'


theta0 and rho0 are basic-state potential temperature and density profiles.

  epsilon_u u - fv + phi_x = 0
  epsilon_v v + fu + phi_y = 0
  epsilon_phi phi + u_x + v_y = -Q_1

the geopotential:

  phi = - (g zT/ pi theta_00) theta
