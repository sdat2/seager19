 # Folder that the ocean model is run from.

```bash
sh ./run-model
```
 
The timings end up getting stored in `timing.txt`.


From what I understand so far, it seems that the 
model parameters are fed in at compilation of the model.


```
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
