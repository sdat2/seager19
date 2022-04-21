#!/bin/bash
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=6666 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5b comp.prwnd=5b
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=SESE comp.sst=5a comp.prwnd=5a
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=SEES  comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=SESS  comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEGG  comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEII  comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEUU  comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,1,2 atm.vary_cloud_const=true,false atm.mem=EEKK  comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=66EE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=666E archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=6EEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a
