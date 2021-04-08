## Folder that the ocean model is run from.

```bash
sh ./run-model.sh
```

The timings end up getting stored in `timing.txt`.

From what I understand so far, it seems that the
model parameters are fed in at compilation of the model.

```bash
 reading file for land/sea mask: DATA/om_mask.nc
 Using a time step of:    10.3655691     hours
 mode=           1   c,dt,hx,hy=   2.98355699      0.307775855      0.307775885      0.307775855    
   uscl, hscl =   2.98355699      0.910139740    
 mode=           2   c,dt,hx,hy=   1.84669042      0.242138833      0.391205281      0.391205251    
   uscl, hscl =   1.84669042      0.348680437    
 saving restart file:  -252.033752    
 saving restart file:  -240.035248    
 saving restart file:  -228.036743    
 saving restart file:  -216.038239    
 saving restart file:  -204.039749    
 saving restart file:  -192.041245    
 saving restart file:  -180.042740    
 saving restart file:  -168.044235    
 saving restart file:  -156.045731    
 saving restart file:  -144.047241    
 saving restart file:  -132.048721    
 saving restart file:  -120.050232    
 saving restart file:  -108.051727    
 saving restart file:  -96.0532227    
 saving restart file:  -84.0547180    
 saving restart file:  -72.0562134    
 saving restart file:  -60.0577087    
 saving restart file:  -48.0592041    
 saving restart file:  -36.0606995    
 saving restart file:  -24.0622101    
 Finished at step =       16903
 <enso time> <  -24.0337982     >
```

Error message for R7
```
*** Error in `../SRC/tcom.exe': realloc(): invalid next size: 0x000000000f921340 ***
======= Backtrace: =========
/lib64/libc.so.6(+0x7f474)[0x7fd9cbe19474]
/lib64/libc.so.6(+0x84861)[0x7fd9cbe1e861]
/lib64/libc.so.6(realloc+0x1d2)[0x7fd9cbe1fe12]
/lib64/libnetcdf.so.7(utf8proc_map+0x89)[0x7fd9ccea6d19]
/lib64/libnetcdf.so.7(utf8proc_NFC+0x23)[0x7fd9ccea6de3]
/lib64/libnetcdf.so.7(NC_findattr+0x31)[0x7fd9cceb8e81]
/lib64/libnetcdf.so.7(+0x6efc1)[0x7fd9cceb8fc1]
/lib64/libnetcdf.so.7(NC3_inq_att+0x24)[0x7fd9cceb91e4]
/lib64/libnetcdf.so.7(nc_inq_att+0x50)[0x7fd9cce9cc30]
/lib64/libnetcdf.so.7(ncattinq+0x2d)[0x7fd9cce9a34d]
/lib64/libnetcdff.so.5(ncainq_+0x7d)[0x7fd9ccbf192d]
../SRC/tcom.exe[0x4310e8]
../SRC/tcom.exe[0x419533]
../SRC/tcom.exe[0x419c88]
../SRC/tcom.exe[0x41f296]
../SRC/tcom.exe[0x404320]
../SRC/tcom.exe[0x4024cd]
/lib64/libc.so.6(__libc_start_main+0xf5)[0x7fd9cbdbc555]
../SRC/tcom.exe[0x40251c]
======= Memory map: ========
00400000-00441000 r-xp 00000000 00:28 37717646892077219                  /home/users/sithom/seager19/ocean/SRC/tcom.exe
00640000-00641000 r--p 00040000 00:28 37717646892077219                  /home/users/sithom/seager19/ocean/SRC/tcom.exe
00641000-00642000 rw-p 00041000 00:28 37717646892077219                  /home/users/sithom/seager19/ocean/SRC/tcom.exe
00642000-0e483000 rw-p 00000000 00:00 0 
0f905000-0f926000 rw-p 00000000 00:00 0                                  [heap]
7fd9bc000000-7fd9bc021000 rw-p 00000000 00:00 0 
7fd9bc021000-7fd9c0000000 ---p 00000000 00:00 0 
7fd9c23bf000-7fd9c6b77000 rw-p 00000000 00:00 0 
7fd9c6b77000-7fd9c6b79000 r-xp 00000000 08:03 1052292                    /usr/lib64/libfreebl3.so
7fd9c6b79000-7fd9c6d78000 ---p 00002000 08:03 1052292                    /usr/lib64/libfreebl3.so
7fd9c6d78000-7fd9c6d79000 r--p 00001000 08:03 1052292                    /usr/lib64/libfreebl3.so
7fd9c6d79000-7fd9c6d7a000 rw-p 00002000 08:03 1052292                    /usr/lib64/libfreebl3.so
7fd9c6d7a000-7fd9c6dda000 r-xp 00000000 08:03 1053022                    /usr/lib64/libpcre.so.1.2.0
7fd9c6dda000-7fd9c6fda000 ---p 00060000 08:03 1053022                    /usr/lib64/libpcre.so.1.2.0
7fd9c6fda000-7fd9c6fdb000 r--p 00060000 08:03 1053022                    /usr/lib64/libpcre.so.1.2.0
7fd9c6fdb000-7fd9c6fdc000 rw-p 00061000 08:03 1053022                    /usr/lib64/libpcre.so.1.2.0
7fd9c6fdc000-7fd9c6fe4000 r-xp 00000000 08:03 1052525                    /usr/lib64/libcrypt-2.17.so
7fd9c6fe4000-7fd9c71e3000 ---p 00008000 08:03 1052525                    /usr/lib64/libcrypt-2.17.so
7fd9c71e3000-7fd9c71e4000 r--p 00007000 08:03 1052525                    /usr/lib64/libcrypt-2.17.so
7fd9c71e4000-7fd9c71e5000 rw-p 00008000 08:03 1052525                    /usr/lib64/libcrypt-2.17.so
7fd9c71e5000-7fd9c7213000 rw-p 00000000 00:00 0 
7fd9c7213000-7fd9c7237000 r-xp 00000000 08:03 1053039                    /usr/lib64/libselinux.so.1
7fd9c7237000-7fd9c7436000 ---p 00024000 08:03 1053039                    /usr/lib64/libselinux.so.1
7fd9c7436000-7fd9c7437000 r--p 00023000 08:03 1053039                    /usr/lib64/libselinux.so.1
7fd9c7437000-7fd9c7438000 rw-p 00024000 08:03 1053039                    /usr/lib64/libselinux.so.1
7fd9c7438000-7fd9c743a000 rw-p 00000000 00:00 0 
7fd9c743a000-7fd9c7456000 r-xp 00000000 08:03 1056270                    /usr/lib64/libsasl2.so.3.0.0
7fd9c7456000-7fd9c7655000 ---p 0001c000 08:03 1056270                    /usr/lib64/libsasl2.so.3.0.0
7fd9c7655000-7fd9c7656000 r--p 0001b000 08:03 1056270                    /usr/lib64/libsasl2.so.3.0.0
7fd9c7656000-7fd9c7657000 rw-p 0001c000 08:03 1056270                    /usr/lib64/libsasl2.so.3.0.0
7fd9c7657000-7fd9c766d000 r-xp 00000000 08:03 1097451                    /usr/lib64/libresolv-2.17.so
7fd9c766d000-7fd9c786d000 ---p 00016000 08:03 1097451                    /usr/lib64/libresolv-2.17.so
7fd9c786d000-7fd9c786e000 r--p 00016000 08:03 1097451                    /usr/lib64/libresolv-2.17.so
7fd9c786e000-7fd9c786f000 rw-p 00017000 08:03 1097451                    /usr/lib64/libresolv-2.17.so
7fd9c786f000-7fd9c7871000 rw-p 00000000 00:00 0 
7fd9c7871000-7fd9c7874000 r-xp 00000000 08:03 1053201                    /usr/lib64/libkeyutils.so.1.5
7fd9c7874000-7fd9c7a73000 ---p 00003000 08:03 1053201                    /usr/lib64/libkeyutils.so.1.5
7fd9c7a73000-7fd9c7a74000 r--p 00002000 08:03 1053201                    /usr/lib64/libkeyutils.so.1.5
7fd9c7a74000-7fd9c7a75000 rw-p 00003000 08:03 1053201                    /usr/lib64/libkeyutils.so.1.5
7fd9c7a75000-7fd9c7a83000 r-xp 00000000 08:03 1097464                    /usr/lib64/libkrb5support.so.0.1
7fd9c7a83000-7fd9c7c83000 ---p 0000e000 08:03 1097464                    /usr/lib64/libkrb5support.so.0.1
7fd9c7c83000-7fd9c7c84000 r--p 0000e000 08:03 1097464                    /usr/lib64/libkrb5support.so.0.1
7fd9c7c84000-7fd9c7c85000 rw-p 0000f000 08:03 1097464                    /usr/lib64/libkrb5support.so.0.1
7fd9c7c85000-7fd9c7c8c000 r-xp 00000000 08:03 1097452                    /usr/lib64/librt-2.17.so
7fd9c7c8c000-7fd9c7e8b000 ---p 00007000 08:03 1097452                    /usr/lib64/librt-2.17.so
7fd9c7e8b000-7fd9c7e8c000 r--p 00006000 08:03 1097452                    /usr/lib64/librt-2.17.so
7fd9c7e8c000-7fd9c7e8d000 rw-p 00007000 08:03 1097452                    /usr/lib64/librt-2.17.so
7fd9c7e8d000-7fd9c80c3000 r-xp 00000000 08:03 1053267                    /usr/lib64/libcrypto.so.1.0.2k
7fd9c80c3000-7fd9c82c3000 ---p 00236000 08:03 1053267                    /usr/lib64/libcrypto.so.1.0.2k
7fd9c82c3000-7fd9c82df000 r--p 00236000 08:03 1053267                    /usr/lib64/libcrypto.so.1.0.2k
7fd9c82df000-7fd9c82ec000 rw-p 00252000 08:03 1053267                    /usr/lib64/libcrypto.so.1.0.2k
7fd9c82ec000-7fd9c82f0000 rw-p 00000000 00:00 0 
7fd9c82f0000-7fd9c8357000 r-xp 00000000 08:03 1063984                    /usr/lib64/libssl.so.1.0.2k
7fd9c8357000-7fd9c8557000 ---p 00067000 08:03 1063984                    /usr/lib64/libssl.so.1.0.2k
7fd9c8557000-7fd9c855b000 r--p 00067000 08:03 1063984                    /usr/lib64/libssl.so.1.0.2k
7fd9c855b000-7fd9c8562000 rw-p 0006b000 08:03 1063984                    /usr/lib64/libssl.so.1.0.2k
7fd9c8562000-7fd9c8569000 r-xp 00000000 08:03 1075878                    /usr/lib64/libaec.so.0.0.10
7fd9c8569000-7fd9c8768000 ---p 00007000 08:03 1075878                    /usr/lib64/libaec.so.0.0.10
7fd9c8768000-7fd9c8769000 r--p 00006000 08:03 1075878                    /usr/lib64/libaec.so.0.0.10
7fd9c8769000-7fd9c876a000 rw-p 00007000 08:03 1075878                    /usr/lib64/libaec.so.0.0.10
Program received signal SIGABRT: Process abort signal.

Backtrace for this error:
#0  0x7FD9CC8D56D7
#1  0x7FD9CC8D5D1E
#2  0x7FD9CBDD044F
#3  0x7FD9CBDD03D7
#4  0x7FD9CBDD1AC7
#5  0x7FD9CBE12F66
#6  0x7FD9CBE19473
#7  0x7FD9CBE1E860
#8  0x7FD9CBE1FE11
#9  0x7FD9CCEA6D18
#10  0x7FD9CCEA6DE2
#11  0x7FD9CCEB8E80
#12  0x7FD9CCEB8FC0
#13  0x7FD9CCEB91E3
#14  0x7FD9CCE9CC2F
#15  0x7FD9CCE9A34C
#16  0x7FD9CCBF192C
#17  0x4310E7 in odb_ifatt_
#18  0x419532 in grids_equiv_
#19  0x419C87 in data_on_model_grid_.part.1 at om_forc.F:?
#20  0x41F295 in qflux_init_
#21  0x40431F in MAIN__ at om_main.F:?
./run-model.sh: line 80: 20041 Aborted                 ../SRC/tcom.exe -i om_run2f -t month.tios
```
