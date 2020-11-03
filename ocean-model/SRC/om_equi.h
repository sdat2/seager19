C	common blocks for EQUIV* routines
C	only exists in EQUIV and routines called from EQUIV
C
      PARAMETER ( MAXZVAL = 1001 )
      
      REAL  NBARCY, NCUR, NBAR
      REAL*8 DZ(MAXZVAL) ,GPINV
      REAL*8  DEQ(MAXZVAL), EV(MAXZVAL), L(MAXZVAL), UEQ(MAXZVAL)
      REAL NSQ(MAXZVAL), Z(MAXZVAL),TEMP1(MAXZVAL)
      INTEGER NOUT
      
C	profile arrays
C	Z	depth
C	DZ	offset delta depth
C	NSQ	buoyancy frequency
C	TEMP1	temperature (potential)
      
      COMMON / BLKEQA / DZ, DEQ, L, UEQ, GPINV
      COMMON / BLKEQR / NBAR, T, GRAV, GRAINV, DML, 
     *                  EDEPTH, PHASEV, RL, RT, FREQSQ, PMODES, Z, TEMP1, NSQ
      COMMON / BLKEQI / NPTS, IPRIME, NMODES , NOUT
c      STATIC EV
      EQUIVALENCE ( EV(1), DEQ(1) )
      
C	GPINV	equivalent gravity for bottom layer
      
