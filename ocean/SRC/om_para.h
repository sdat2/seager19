C  meqsf.h  parameters for multimode model
C  here are four parameters which are almost mnemonic.

C  MAXMODES=MMODES is the maximum no. of modes.
      INTEGER MAXMODES, NXMAX, NYMAX, MAXDIM, MAXNZ, NZMAX
      INTEGER MXY, MMODES
      PARAMETER ( MAXMODES= 20 )

C  These two parameters define the maximum grid size in the ocean basin.
      PARAMETER ( NXMAX=400 )
      PARAMETER ( NYMAX=320 )
C	a parameter for the greater of the two (NYMAX or NXMAX)
	PARAMETER (MAXDIM=NYMAX)

      PARAMETER ( MAXNZ  = 51)

      PARAMETER ( NZMAX  = MAXNZ )
      PARAMETER ( MMODES = MAXMODES )
      PARAMETER ( MXY  = NXMAX*NYMAX )

