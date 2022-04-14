#!/bin/bash

python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=E6E6 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5c comp.prwnd=5c
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=E666 archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5b comp.prwnd=5b
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=E66E archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a
python src/main.py -m atm.e_frac=0.5,2 atm.vary_cloud_const=true,false atm.mem=EEEE archive_dir=/gws/nopw/j04/ai4er/users/sdat2/rep comp.sst=5a comp.prwnd=5a
