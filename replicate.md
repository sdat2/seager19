# Replicate each figure

```bash


python src/main.py name=std_uncoup coup.iterations=1 atm.mem=EEEf atm.vary_cloud_const=true archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false

python src/main.py name=std_coup coup.iterations=6 atm.mem=EEEf atm.vary_cloud_const=true archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false


python src/main.py name=CMIP5_coup coup.iterations=6 atm.mem=EECC atm.vary_cloud_const=true archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false


python src/main.py name=CR_coup coup.iterations=6 atm.mem=EEEC atm.vary_cloud_const=true archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep wandb=false

# The effective relative humidity, r, is assumed to be spatially 
# uniform at 0.80 in our standard model
# In considering the cause of bias in CMIP5 models, we instead imposed 
# the spatially varying climatological annual mean r from, first, 
# ECMWF and, second, the CMIP5 multimodel mean.
# EEEc, EEEf
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
# c - forced with atmosphere-ocean
# 5 - coupled models
# a - CM-ECMWF world
# b - CM-ECMWF C-RH - change the relative humidity to CMIP5
# c - CM-ECMWF C-RH W - change the relative humidity and the wind to CMIP5
# d - CM-CMIP5 world qfluxed towards the CMIP5 mmm
# f - unknown thermocline thing.
# g - ECMWF graph
# h - CMIP5 multimodel mean
```
