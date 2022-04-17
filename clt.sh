#!/bin/bash

python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=ESES comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=ESSS comp.sst=5b comp.prwnd=5b
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=ESSE comp.sst=5a comp.prwnd=5a
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=EEEE comp.sst=5a comp.prwnd=5a
