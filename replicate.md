# Replicate each figure

```bash

python src/main.py name=test25

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true atm.mem=EESS comp.sst=5c comp.prwnd=5c


python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false coup.stress_trend=false,true coup.add_stress=false,true atm.prcp_land=1,0 ocean.flux_once=true atm.mem=EEEE  comp.sst=5a comp.prwnd=5a oc.f1prime=-0.0021

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false coup.stress_trend=false,true coup.add_stress=false,true atm.prcp_land=1,0 ocean.flux_once=true atm.mem=EEEE  comp.sst=5a comp.prwnd=5a oc.f1prime=-0.006

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false coup.stress_trend=false,true coup.add_stress=false,true atm.prcp_land=1,0 ocean.flux_once=true atm.mem=EEEE  comp.sst=5a comp.prwnd=5a oc.f1prime=-0.003

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false coup.stress_trend=false,true coup.add_stress=false,true atm.prcp_land=1,0 ocean.flux_once=true atm.mem=EEEE  comp.sst=5a comp.prwnd=5a oc.f1prime=-0.004

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false coup.stress_trend=false,true coup.add_stress=false,true atm.prcp_land=0 ocean.flux_once=true,false atm.mem=EEEE  comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true atm.mem=EEEE comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EECCE comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEEEE comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEGG comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEUU comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEKK comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEII comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=false,true atm.mem=SSSS comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEEC comp.sst=5b comp.prwnd=5b

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=false,true atm.mem=EEES comp.sst=5b comp.prwnd=5b

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EESE comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EESS comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=false,true atm.mem=SSSSSSE comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=false,true atm.mem=SSSSSSS comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=false,true atm.mem=EEES comp.sst=5b comp.prwnd=5b

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false  atm.mem=EESE comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=false,true atm.mem=EEEC comp.sst=5b comp.prwnd=5b

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false  atm.mem=EECE comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false coup.add_stress=true coup.stress_trend=true  atm.mem=EEEE comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false coup.add_stress=true coup.stress_trend=true  atm.mem=EECC comp.sst=5c comp.prwnd=5c

python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true atm.mem=EE6E archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=false atm.mem=EE6E archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a

python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false  atm.mem=EEEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a

# explore cd sensitivity:

python src/main.py -m coup.c_d=1.25e-3,1.5e-3,1.75e-3,2e-3,2.25e-3,2.5e-3,2.75e-3,3e-3 atm.mem=EECC atm.e_frac=2 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5c comp.prwnd=5c

python src/main.py -m coup.c_d=1.25e-3,1.5e-3,1.75e-3,2e-3,2.25e-3,2.5e-3,2.75e-3,3e-3 atm.mem=EEEC atm.e_frac=2 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5b comp.prwnd=5b

python src/main.py -m coup.c_d=1.25e-3,1.5e-3,1.75e-3,2e-3,2.25e-3,2.5e-3,2.75e-3,3e-3 atm.mem=EECE atm.e_frac=2 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5a comp.prwnd=5a

python src/main.py -m coup.c_d=1.25e-3,1.5e-3,1.75e-3,2e-3,2.25e-3,2.5e-3,2.75e-3,3e-3 atm.mem=EEEE atm.e_frac=2 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5a comp.prwnd=5a

## add more runs.

python src/main.py -m coup.c_d=0.5e-3,3.25e-3 atm.mem=EECC atm.e_frac=0.5 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5c comp.prwnd=5c

python src/main.py -m coup.c_d=0.5e-3,0.75e-3,2.75e-3,3.25e-3  atm.mem=EEEC atm.e_frac=0.5 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5b comp.prwnd=5b

python src/main.py -m coup.c_d=0.5e-3,2.75e-3,3.25e-3  atm.mem=EECE atm.e_frac=2 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5a comp.prwnd=5a

python src/main.py -m coup.c_d=0.5e-3,0.75e-3,2.75e-3,3.25e-3  atm.mem=EEEE atm.e_frac=2 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/all comp.sst=5a comp.prwnd=5a

## replication of parameters in paper

python src/main.py name=N_std_uncoup coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d

python src/main.py name=N_uncoup_noheat coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false atm.prcp_land=0 comp.htherm=4b comp.sst=1d comp.prwnd=2c

python src/main.py name=N_ECMWF_uncoup atm.mem=EEEE coup.iterations=1 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d

python src/main.py name=N_std_coup atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep

python src/main.py name=N_ECMWF_coup atm.mem=EEEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a

python src/main.py name=N_C_RH_W2_coup atm.mem=EECC archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5c comp.prwnd=5c

python src/main.py name=N_C_RH_coup atm.mem=EEEC  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5b comp.prwnd=5b

python src/main.py name=N_C_W_coup atm.mem=EECE  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a


# making the cloud constant constant  atm.vary_cloud_const=false

python src/main.py name=A_std_uncoup coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d atm.vary_cloud_const=false

python src/main.py name=A_uncoup_noheat coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false atm.prcp_land=0 comp.htherm=4b comp.sst=1d comp.prwnd=2c atm.vary_cloud_const=false

python src/main.py name=A_ECMWF_uncoup atm.mem=EEEE coup.iterations=1 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d atm.vary_cloud_const=false

python src/main.py name=A_std_coup atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep atm.vary_cloud_const=false

python src/main.py name=A_ECMWF_coup atm.mem=EEEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a atm.vary_cloud_const=false

python src/main.py name=A_C_RH_W2_coup atm.mem=EECC archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5c comp.prwnd=5c atm.vary_cloud_const=false

python src/main.py name=A_C_RH_coup atm.mem=EEEC  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5b comp.prwnd=5b atm.vary_cloud_const=false

python src/main.py name=A_C_W2_coup atm.mem=EECE  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a atm.vary_cloud_const=false
 
# Change the epsilon fraction.

python src/main.py name=E2_std_uncoup coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d atm.e_frac=2

python src/main.py name=E2_uncoup_noheat coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false atm.prcp_land=0 comp.htherm=4b comp.sst=1d comp.prwnd=2c atm.e_frac=2

python src/main.py name=E2_ECMWF_uncoup atm.mem=EEEE coup.iterations=1 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d atm.e_frac=2

python src/main.py name=E2_std_coup atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep atm.e_frac=2

python src/main.py name=E2_ECMWF_coup atm.mem=EEEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a atm.e_frac=2

python src/main.py name=E2_C_RH_W_coup atm.mem=EECC archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5c comp.prwnd=5c atm.e_frac=2

python src/main.py name=E2_C_RH_coup atm.mem=EEEC  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5b comp.prwnd=5b atm.e_frac=2

python src/main.py name=E2_C_W_coup atm.mem=EECE  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a atm.e_frac=2


# Change both.

python src/main.py name=AE_std_uncoup coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d atm.e_frac=2 atm.vary_cloud_const=false

python src/main.py name=AE_uncoup_noheat coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false atm.prcp_land=0 comp.htherm=4b comp.sst=1d comp.prwnd=2c atm.e_frac=2 atm.vary_cloud_const=false

python src/main.py name=AE_ECMWF_uncoup atm.mem=EEEE coup.iterations=1 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d atm.e_frac=2 atm.vary_cloud_const=false

python src/main.py name=AE_std_coup atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep atm.e_frac=2 atm.vary_cloud_const=false

python src/main.py name=AE_ECMWF_coup atm.mem=EEEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a atm.e_frac=2 atm.vary_cloud_const=false

python src/main.py name=AE_C_RH_W_coup atm.mem=EECC archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5c comp.prwnd=5c atm.e_frac=2 atm.vary_cloud_const=false

python src/main.py name=AE_C_RH_coup atm.mem=EEEC  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5b comp.prwnd=5b atm.e_frac=2 atm.vary_cloud_const=false

python src/main.py name=AE_C_W_coup atm.mem=EECE  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a atm.e_frac=2 atm.vary_cloud_const=false

# The effective relative humidity, r, is assumed to be spatially 
# uniform at 0.80 in our standard model
# In considering the cause of bias in CMIP5 models, we instead imposed 
# the spatially varying climatological annual mean r from, first, 
# ECMWF and, second, the CMIP5 multimodel mean.
# EEEC, EEEf 
# The 'standard model' from section 3 uses
# 1 - forced ocean
# d - rising CO2, observed winds
# e - rising CO2, fixed winds
# f - fixed CO2, rising winds
# 2 - forced atmosphere
# c - no heating over land, ECMWF forcing
# d - heating over land, ECMWF forcing
# 3 - coupled atmosphere ocean model
#     heating over land, ECMWF inputs
# a - sst change
# b - prcp, utrend, vtrend change
# 4 - Trend in the thermocline
# a - ORAS4 model
# b - forced with ORAS4 winds
# c - coupled atmosphere-ocean
# 5 - coupled models
# a - CM-ECMWF world
# b - CM-ECMWF C-RH - change the relative humidity to CMIP5
# c - CM-ECMWF C-RH W - change the relative humidity and the wind to CMIP5
# d - CM-CMIP5 world qfluxed towards the CMIP5 mmm
# f - unknown thermocline thing.
# g - ECMWF graph
# h - CMIP5 multimodel mean
```
