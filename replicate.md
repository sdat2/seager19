# Replicate each figure

```bash


python src/main.py name=std_uncoup coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d # set off

python src/main.py name=uncoup_noheat coup.iterations=1 atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false atm.prcp_land=0 comp.htherm=4b comp.sst=1d comp.prwnd=2c # set off

python src/main.py name=ECMWF_uncoup atm.mem=EEEE coup.iterations=1 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false comp.htherm=4b comp.sst=1d comp.prwnd=2d

python src/main.py name=std_coup atm.mem=EEEf archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep

python src/main.py name=ECMWF_coup atm.mem=EEEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a

python src/main.py name=CRH_W_coup atm.mem=EECC archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5c comp.prwnd=5c

python src/main.py name=CRH_coup atm.mem=EEEC  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5b comp.prwnd=5b

python src/main.py name=CW_coup atm.mem=EECE  archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a

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
