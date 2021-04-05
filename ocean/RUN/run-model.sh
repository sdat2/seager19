# usage: sh ./run-model.sh

# STEP 1 -  spinup the dynamics - ignore the sst equation and use climo windstress

echo "A: $(date)" >> timing.txt

../SRC/tcom.exe -i om_spin -t spin.tios  # 6 minutes

echo "B: $(date)" >> timing.txt

../SRC/tios2cdf.exe -f output/om_spin # 3 seconds

echo "C: $(date)" >> timing.txt

rm -rf output/om_spin.data output/om_spin.indx    # < 1 seconds
cp -f output/om_spin.save output/om_spin.20y.restart   # < 1 seconds

echo "D: $(date)" >> timing.txt  

# STEP 2 -  diagnose the qflux needed to reproduce sst climatology

../SRC/tcom.exe -i om_diag -t diag.tios # 55 seconds

echo "E: $(date)" >> timing.txt  

../SRC/tios2cdf.exe -f output/om_diag # < 1 seconds

echo "F: $(date)" >> timing.txt

rm -rf output/om_diag.data output/om_diag.indx # < 1 seconds
cp -f output/om_diag.save output/om_diag.2y.restart # < 1 seconds

# if you have ingrid installed, do the next line,
# otherwise take the last year of qfluxes from om_diag.nc and 
# pop it into DATA/qflx.nc with your favorite tool
# /usr/local/bin/ingrid qflx.ing


# STEP 3 -  run the dynamics + sst in full mode

echo "G: $(date)" >> timing.txt

../SRC/tcom.exe -i om_run2f -t month.tios    # 25 minutes 

echo "H: $(date)" >> timing.txt

../SRC/tios2cdf.exe -f output/om_run2f     # 9 seconds

echo "I: $(date)" >> timing.txt

rm -rf output/om_run2f.data output/om_run2f.indx     # < 1 seconds

echo "J: $(date)" >> timing.txt
