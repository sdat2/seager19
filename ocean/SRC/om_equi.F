*----------------------------------------------------------------------
      subroutine VMODESSUB(NIN,PIN,ANBV2,TEMP,NUMMOD,ZCUT,ZIS,ZTOP,VGRID,
     *     NZ,hstarp,CPH,LIGHT,SLVL,HVM,VVM,HCM,HVMSF,ahmix0,ahmixp)
      
      INCLUDE 'om_equi.h'
*----------------------------------------------------------------------
C NAME
C	vmodessub(9)		- subroutine to calculate vertical modes
C SYNOPSIS
C	subroutine VMODESSUB(NIN,PIN,ANBV2,TEMP,NUMMOD,ZCUT,
C     *  ZIS,ZTOP,VGRID,NZ,hstarp,CPH,LIGHT,SLVL,HVM,VVM,HCM)
C Inputs
C	NIN	[]	number of points in input profile
C	PIN	[P>	input pressures (or depths (m))
C	ANBV2	[ANBV2>	input BV frequency squared (rad/s)squared
C	TEMP	[TEMP>	input temperature profile
C	NUMMOD	()	number of modes to calculate
C	ZCUT		equivalent bottom layer depth parameter
C			has no effect if greater than water col depth
C	ZIS	char	vertical coordinate name
C			either Depth or Temperature
C	ZTOP		smallest coordinate to interpolate output to
C	VGRID		output grid
C	NZ	{}	number of points in output grid
C Outputs
C	hstarp   	1./(total Lighthill coefficient over NUMMOD modes)
C	CPH	(CPH>	phase speeds (m/s)
C	LIGHT	(LIGH>	Lighthill coefficents (1/m)
C	SLVL	(SLVL>	Sea Level (m)
C	HVM	{HVM)	Horizontal Velocity Modes
C	VVM	{VVM)	Vertical Velocity Modes
C	HCM	{HCM)	Heat Content Modes
C       HVMSF           surface value of HVM
C	ahmix0		depth of 20degree isotherm
C	ahmixp		value of vertical velocity mode at 20degree isotherm
C
C	all arguments are real*4
C	double precision imsl library required
C
C DESCRIPTION
C	This program uses a tridiagonal matrix solver to calculate dynamical
C	modes from a buoyancy frequency profile.  The problem is solved with
C	a free surface boundary condition.
C BUGS
C	Clearly this is an awkward subroutine, with way too many arguments.
C AUTHOR
C       Benno provided this version of the equiv* subroutines
*----------------------------------------------------------------------

      REAL	PIN(NIN),ANBV2(NIN),TEMP(NIN)
      REAL	CPH(NUMMOD),LIGHT(NUMMOD),SLVL(NUMMOD)
      REAL	HVM(NZ,NUMMOD),VVM(NZ,NUMMOD),HCM(NZ,NUMMOD)
      REAL	ahmix0(NUMMOD),ahmixp(NUMMOD),hvmsf(nummod)
      
      REAL	ZCUT,ZTOP,VGRID(NZ),hstarp
      
      CHARACTER*(*) ZIS
      CHARACTER*256	FNAME
      CHARACTER*80	LINE
      
C	altered 27 April 1987 to include from local directory
C	19 May 1987 various modifications including equivalent bottom layer,
C		double precision eigenvalue/vector calculation, temperature
C		coordinates, and new printout
C	altered 15 AUG 1997 - Naomi- arbitrary output grid, vgrid

      GRAV   = 9.78049
      GRAINV = 1.0 / GRAV
      
C DML is the mixed layer depth. disabled!!!!!!!!!!!!
C      DML = HSFC
C 	NOUT is the COMMON block copy of NZ
      NOUT = NZ
C NMODES includes the barotropic mode (#1) which will be discarded at the end.
C NUMMOD includes only baroclinic modes, (#2 to #n).
      NMODES = NUMMOD + 1
      
C  Read the water profile and...
      CALL EQUIV1(ZCUT,ZIS,ZTOP,VGRID,NZ,NIN,PIN,ANBV2,TEMP)
      
C  Compute the tridiagonal matrix
      CALL EQUIV2
      
C  Call the eigenvalue routine
      CALL EQUIV3
      
C  calculate eigenvectors, fns of eigenvalues
      CALL EQUIV4(ZIS,ZTOP,VGRID,NZ,hstarp,
     *     CPH,LIGHT,SLVL,HVM,VVM,HCM,hvmsf,ahmix0,ahmixp)
C	HVM	{HVM)	selected horizontal velocity mode values
C	VVM	{VVM)	selected vertical velocity mode values
C	HCM	{HCM)	selected heat content mode values
      
      return
      END
      
*----------------------------------------------------------------------
      subroutine MODEPA(nx,ny,nmodes,npt,time_mul,delt,delx,dely,xm,ym,
     *            X_MIN,X_MAX,Y_MIN,Y_MAX,cph,phasev,cl,sealvl,dtm,A_rayl,rayl,
     *            A_lapl,rlap,hxm,hym,fcor,fc,hscl,uscl,dscl,ahmixp,
     *            scale_L)
      implicit none
*----------------------------------------------------------------------
C manpage modepa
C NAME
C	modepa(9)[vmodessub]	- subroutine to calculate scales for modes
C SYNOPSIS
C	subroutine MODEPA([grid_stuff],PHASEV,SEALVL,CL,HSCL,USCL,DSCL)
C Inputs
C       nx              number of zonal grid points 
C       ny              number of meridional grid points 
C       nmodes          number of vertical modes
C       time_mul        1: reset delt = time_mul*delx 
C       delt            time step (s)
C       A_rayl          Rayleigh friction coefficient
C       A_lapl          Laplacian diffusion coefficient
C       X_MIN,X_MAX     zonal limits (degrees)
C       Y_MIN,Y_MAX     meridional limits (degrees)
C	phasev		phase velocity (m/s)
C	sealvl		sea level coefficient
C	cl 		Lighthill coupling coef.:	du/dt + ... = taux/CL
C	ahmixp 		vertical velocity eigenfunction at z(temp=20)
C Outputs
C       delt            modified time step (s)
c       delx            zonal grid spacing (m)
c       dely            merid grid spacing (m)
c       xm              zonal grid (m)
c       ym              merid grid (m)
C       dtm             scaled time
C       rayl            scaled rayleigh friction
C       rlap            scaled laplacian coefficient
C       hxm             scaled zonal grid spacing
C       hym             scaled merid grid spacing
C       fc              scaled coriolis parameter
C       fcor            unscaled coriolis parameter (1/s)
C	hscl		scaling for pressure (m)
C	uscl		scaling for velocity (m/s)
C	dscl		scaling for vertical displacement (m)
C	ahmixp 		scaled ahmixp for tdeep
C DESCRIPTION
C	modepa is a subroutine to calculate mode scales from the information
C	that vmodessub returns.
C BUGS
C	It has not been tested yet
C AUTHOR
C       Benno provided the skeleton of this subroutine
*----------------------------------------------------------------------

      integer nx,ny,nmodes, npt
      real delt,delx,dely,time_mul,X_MIN,X_MAX,Y_MIN,Y_MAX,
     *     A_rayl,A_lapl
      real phasev(nmodes),sealvl(nmodes,nx*ny),cl(nmodes,nx*ny)
     *         ,dtm(nmodes),hxm(nmodes),hym(nmodes)
     *         ,dscl(nmodes,nx*ny),uscl(nmodes),hscl(nmodes)
     *         ,fc(ny,nmodes),ahmixp(nmodes,nx*ny),xm(nx),ym(ny),fcor(ny)
     *         ,rayl(nmodes),rlap(nmodes),cph(nmodes,nx*ny)

      real grav, grainv, r_earth, conv, rt, rl, to_rad, CNST_2OMEGA, y_ave
      real phi_0, CNST_BETA, CNST_F0, yconv, xconv, xmin, ymin 
      real scale_L
      integer mode, i, j, ipt

      logical use_profile
      integer nvert_prof
      real xlon_prof, ylat_prof, hmix_param
      integer ipt_prof
      common /vert_profile/ use_profile,xlon_prof, ylat_prof, nvert_prof,
     *                      ipt_prof, hmix_param
      common /constants/ xconv,yconv,rl

      grav   = 9.780490
      grainv = 1.0 / grav
      r_earth = 6379000.
      to_rad = atan(1.)/45.
      conv = r_earth* to_rad

      CNST_2OMEGA = 2.* 2.* (4.*atan(1.))/ (24.* 3600.)! 2*2*pi/day in 1/s
      y_ave = (Y_MIN + Y_MAX)/2.0
      PHI_0 = to_rad*y_ave

      CNST_BETA   = CNST_2OMEGA * cos(PHI_0)/R_EARTH
      CNST_F0     = CNST_2OMEGA * ( sin(PHI_0) - PHI_0 * cos(PHI_0) )

      yconv = conv
      xconv = conv*cos(PHI_0)

      delx  = xconv*(X_MAX - X_MIN)/(nx - 1)   !  x spacing in meters
      dely  = yconv*(Y_MAX - Y_MIN)/(ny - 1)   !  y spacing in meters
      xmin  = xconv*X_MIN                      !  first x value (m)
      ymin  = yconv*Y_MIN                      !  first y value (m)

      do i = 1, nx
         xm(i) = xmin + (i-1)*delx
      enddo
      do j = 1, ny
         ym(j) = ymin + (j-1)*dely
      enddo

      do j = 1, ny
         fcor(j) = CNST_F0 + CNST_BETA*ym(j)
      enddo

      do mode = 1, nmodes
         phasev(mode) = cph(mode,ipt_prof)
      enddo

      if (time_mul.ne.0.) delt = time_mul*delx/phasev(1)

      print*,'Using a time step of: ',delt/3600,'hours'

      do mode = 1, nmodes
         rl = sqrt( phasev(mode) / CNST_BETA )
         rt = sqrt( 1.0 / ( phasev(mode) * CNST_BETA ) )

C               Rayleigh friction, scaled for each mode
         rayl(mode) = rt * A_rayl * (phasev(1)/phasev(mode))**2
         rlap(mode) = rt * A_lapl / rl**2

C       	scale for displacement (isopycnal height) (m)
         do i = 1, nx*ny
            dscl(mode,i)= grainv* rl/ (cl(mode,i)* sealvl(mode,i) )
            
c precompute coefs for tdeep:
            ahmixp(mode,i) =-dscl(mode,i) * ahmixp(mode,i) 
         enddo
            
C	        scale for velocity  (u,v) (m/s)
         uscl(mode) = phasev(mode)
      
C	        scale for pressure perturbations (h) (m)
         hscl(mode) = grainv * (phasev(mode)**2)

C	        scaled time step
         dtm(mode) = delt/rt

C	        scaled zonal grid spacing
         hxm(mode) = delx/rl

C	        scaled meridional grid spacing
         hym(mode) = dely/rl

         print*,'mode=',mode,'  c,dt,hx,hy=',phasev(mode),dtm(mode),hxm(mode),hym(mode)
         print*,'  uscl, hscl =',uscl(mode),hscl(mode)

         do j = 1, ny
            fc(j,mode) = ym(j)/rl + CNST_F0*rt
         enddo

      enddo

      scale_L = sqrt( phasev(1) / CNST_BETA )

      return
      end
      
*----------------------------------------------------------------------
      subroutine single_mode(nx,ny,time_mul,delt,delx,dely,xm,ym,
     *          X_MIN,X_MAX,Y_MIN,Y_MAX,phasev,H,H1,hvm,hvmsf,dtm,A_rayl,rayl
     *          ,hxm,hym,fcor,fc,dscl,hscl,uscl,ahmix0,ahmixp,hstarp,area,
     *          scale_L)
      implicit none
*----------------------------------------------------------------------
C Inputs
C       nx              number of zonal grid points 
C       ny              number of meridional grid points 
C       time_mul        1: reset delt = time_mul*delx 
C       delt            time step (s)
C       A_rayl          Rayleigh friction coefficient
C       X_MIN,X_MAX     zonal limits (degrees)
C       Y_MIN,Y_MAX     meridional limits (degrees)
C	phasev		phase velocity (m/s)
C	H		depth of active ocean
C	H1		depth of surface layer
C Outputs
C       delt            modified time step (s)
c       delx            zonal grid spacing (m)
c       dely            merid grid spacing (m)
c       xm              zonal grid (m)
c       ym              merid grid (m)
C       dtm             scaled time
C       rayl            scaled rayleigh friction
C       hxm             scaled zonal grid spacing
C       hym             scaled merid grid spacing
C       fc              scaled coriolis parameter
C       fcor            unscaled coriolis parameter (1/s)
C	hscl		scaling for pressure (m)
C	uscl		scaling for velocity (m/s)
C	hvm		eigenfunction reset to 1.
C	hvmsf		surface value eigenfunction reset to 1.
C	ahmix0 		rescaled ahmix0
C	ahmixp 		rescaled ahmixp for tdeep
*----------------------------------------------------------------------

      integer nx,ny
      real delt,delx,dely,time_mul,X_MIN,X_MAX,Y_MIN,Y_MAX,A_rayl, area
      real phasev,dtm,hxm,hym, h, h1, hvm(1), hvmsf(1)
     *         ,dscl(1),uscl,hscl,hstarp(1)
     *         ,fc(ny),ahmix0(1),ahmixp(1),xm(nx),ym(ny),fcor(ny),rayl(1)

      real grav, grainv, r_earth, conv, rt, rl, to_rad, CNST_2OMEGA, y_ave
      real phi_0, CNST_BETA, CNST_F0, yconv, xconv, xmin, ymin, scale_L
      integer mode, i, j, ipt
      common /constants/ xconv,yconv,rl

      grav   = 9.780490
      grainv = 1.0 / grav
      r_earth = 6379000.
      to_rad = atan(1.)/45.
      conv = r_earth* to_rad

      CNST_2OMEGA = 2.* 2.* (4.*atan(1.))/ (24.* 3600.)! 2*2*pi/day in 1/s
      y_ave = (Y_MIN + Y_MAX)/2.0
      PHI_0 = to_rad*y_ave

      CNST_BETA   = CNST_2OMEGA * cos(PHI_0)/R_EARTH
      CNST_F0     = CNST_2OMEGA * ( sin(PHI_0) - PHI_0 * cos(PHI_0) )

      yconv = conv
      xconv = conv*cos(PHI_0)

      delx  = xconv*(X_MAX - X_MIN)/(nx - 1)   !  x spacing in meters
      dely  = yconv*(Y_MAX - Y_MIN)/(ny - 1)   !  y spacing in meters
      xmin  = xconv*X_MIN                      !  first x value (m)
      ymin  = yconv*Y_MIN                      !  first y value (m)

      do i = 1, nx
         xm(i) = xmin + (i-1)*delx
      enddo
      do j = 1, ny
         ym(j) = ymin + (j-1)*dely
      enddo

      do j = 1, ny
         fcor(j) = CNST_F0 + CNST_BETA*ym(j)
      enddo

      if (time_mul.ne.0.) delt = time_mul*delx/phasev

      rl = sqrt( phasev / CNST_BETA )
      scale_L = rl
      rt = sqrt( 1.0 / ( phasev * CNST_BETA ) )
      
      rayl(1) = rt * A_rayl
      
C	        scale for velocity  (u,v) (m/s)
      uscl = phasev
      
C	        scale for pressure perturbations (h) (m)
      hscl = H
      hvm(1) = 1.
      hvmsf(1)= 1.

      ahmix0(1) = H
      ahmixp(1) = hscl
      
C	        scaled time step
      dtm = delt/rt
      
C	        scaled zonal grid spacing
      hxm = delx/rl
      
C	        scaled meridional grid spacing
      hym = dely/rl

      area = rt/hscl/delx/dely
      
c         print*,'mode',mode,phasev(mode),dtm(mode),hxm(mode),hym(mode)

      print*,'rayl,hscl,rl,rt: ',rayl(1),hscl,rl,rt
      
      do j = 1, ny
         fc(j) = (j-1)*hym + ymin/rl + CNST_F0*rt
      enddo

      dscl(1) = 1.
      do ipt = 2, nx*ny
         dscl(ipt)  = dscl(1)
         hvmsf(ipt) = hvmsf(1)
         ahmixp(ipt)= ahmixp(1)
         ahmix0(ipt) = ahmix0(1)
         hvm(ipt)    = hvm(1)
         hstarp(ipt) = hstarp(1)
      enddo
      
      return
      end
      
*----------------------------------------------------------------------

      subroutine EQUIV1(ZCUT,ZIS,ZTOP,VGRID,NZ,NPTSI,PIN,ANBV2,TEMP)
      
      INCLUDE 'om_equi.h'
*----------------------------------------------------------------------
C	equiv1[vmodes]	- profile processor for vmodes
C Inputs
C	NPTS	[]	number of points in input profile
C	PIN	[P>	input pressures
C	ANBV2	[ANBV2>	input BV frequency squared (rad/s)squared
C	TEMP	[TEMP>	input temperature profile
C Outputs
C	NBAR	common	Vertically averaged buoyancy frequency in rad/sec
C	Z	common  NPTS array of depths
C	NSQ	common  NPTS array of squared buoyancy frequency
      
      CHARACTER ZIS*(*)
      REAL PIN(NPTSI),ANBV2(NPTSI),TEMP(NPTSI), VGRID(NZ)
      
C	Internal arrays
C	EDENS1	equivalent density
C	HC0	basic state heat content
      
      REAL EDENS1(MAXZVAL),HC0(MAXZVAL)
      
      INTEGER	HEQLUN
C	variables to pass grids to LININT
      REAL VIN(MAXZVAL)
      
      NPTS = NPTSI

      DO 99 J = 1 , NPTS
         Z(J) = PIN(J)
         NSQ(J) = AMAX1(0.,ANBV2(J))
         TEMP1(J) = TEMP(J)
   99 CONTINUE
      
C     *  checks whether first point in profile is at Z=0   *
      IF ( Z(1) .LT. 0.01 ) GO TO 350
      NPTS = NPTS+1
      I1   = NPTS
      NPTSI = NPTS
      
      DO 330 I = 2, NPTS
         Z(I1)  =   Z(I1-1)
         NSQ(I1)  = NSQ(I1-1)
         TEMP1(I1) = TEMP1(I1-1)
         I1   = I1-1
  330 CONTINUE
      Z(1) = 0.0
      
  350 CONTINUE
      
C	computes depth differences
      DO 410 I = 2, NPTS
         DZ(I) = Z(I) - Z(I-1)
  410 CONTINUE
      DZ(1) = Z(1)
      
      IMAX = NPTS - 1
      NBAR = SQRT(NSQ(1))*Z(2) + SQRT(NSQ(NPTS))*(Z(NPTS)-Z(NPTS-1))
      DO 420 I = 2, IMAX
         NBAR = NBAR + SQRT(NSQ(I))*(Z(I+1)-Z(I-1))
  420 CONTINUE
      NBAR = NBAR / (2*Z(NPTS))
      
C	computes equivalent density (integrated buoyancy frequency)
      SUM = 0.
      EDENS1(1) = 1.
      DO 500 I = 2,NPTS
         SUM = SUM + (NSQ(I)+NSQ(I-1))*DZ(I)*.5*GRAINV
         EDENS1(I) = EXP(SUM)
  500 CONTINUE
      
C	computes basic state heat content
      SUM = 0.
      HC0(1) = SUM
      DO 505 I = 2 , NPTS
         SUM = SUM + (TEMP1(I)+TEMP1(I-1))*DZ(I)*.5
         HC0(I) = SUM
  505 CONTINUE
      
C	uses full profile for TEMP and EDENS1 interpolation
C	even if profile gets cut
      
C	decides on grid for interpolated variables
      IF(ZIS(1:1).EQ.'T')THEN
C	temperature grid
         DO 510 IZ = 1 , NPTS
            VIN(IZ) = TEMP1(IZ)
  510    CONTINUE
      ELSE
C	depth grid
         DO 520 IZ = 1 , NPTS
            VIN(IZ) = Z(IZ)
  520    CONTINUE
      ENDIF
C	output grid is based on ZTOP and ZBOT
      ZVAL = ZTOP
      
C-----------------------------------------------------------------------
C	modifies profile with an equivalent bottom layer if requested
C	if ZCUT < Z(NPTS) then creates layer
      
      IF(ZCUT.GT.Z(NPTS))THEN
         GPINV = 0.0
         DZ(NPTS+1) = 1
      ELSE
C	rounds ZCUT to nearest gridpoint
         CALL BINSERCH ( ZCUT, Z, NPTS, IPT )
         IF(  (ZCUT-Z(IPT)) .GT. (Z(IPT+1)-ZCUT) )IPT = IPT + 1
         ZCUT = Z(IPT)
         
C	 computes inverse reduced gravity
         SUM = 0.0
         IMIN = IPT + 1
         IL = IPT
         DO 600 I = IMIN , NPTS
            SUM = SUM + (NSQ(I) + NSQ(IL))*DZ(I)
            IL = I
  600    CONTINUE
         GPINV = 2.0 / SUM
         
C	patchs bottom of profile
         DZ(IPT+1) 	= Z(NPTS) - ZCUT
         Z(IPT+1) 	= Z(NPTS)
         NSQ(IPT+1)	= NSQ(NPTS)
         TEMP1(IPT+1) 	= TEMP1(NPTS)
         
C	remaining patches unecessary since arrays are local to EQUIV1
         EDENS1(IPT+1)	= EDENS1(NPTS)
         HC0(IPT+1)	= HC0(NPTS)
         
         NPTS = IPT + 1
      END IF

      return
      END
      
*----------------------------------------------------------------------
      subroutine BINSERCH ( ZVAL, Z, NPTS, IPT )
*----------------------------------------------------------------------
C
C  03 June 1986
C  Search the Z array and get an IPT value.
C  We assume the Z array is sorted in ascending order.
C  If ZVAL is not between Z(1) and Z(NPTS) then IPT = 0.
C


      DIMENSION Z(NPTS)
      
C  In range?
      IF ( ZVAL .LT. Z(1) .OR. ZVAL .GT. Z(NPTS) ) THEN
         IPT = 0
         return
      END IF
      
C Start and end pointers, which will be 'squeezed together'.
      IS = 1
      IE = NPTS
      
C A 'WHILE' statement would be very nice at this point.
   40 CONTINUE
      
C The pointer and the spread/difference between the bounds
      IPT = ( IS + IE ) / 2
      IDIF = IE - IS
      
C A direct hit!
      IF ( ZVAL .EQ. Z(IPT) ) return
      
C Our bounds can't be squeezed any more.
      IF ( IDIF .EQ. 1 ) THEN
         IPT = IS
         return
      END IF
      
C Determine which bound to 'move in'.
      IF ( ZVAL .LT. Z(IPT) ) THEN
         IE = IPT
      ELSE IF ( ZVAL .GT. Z(IPT) ) THEN
         IS = IPT
      END IF
      
C Do it again.
      GO TO 40
      
      END

*----------------------------------------------------------------------
      subroutine EQUIV2
      INCLUDE 'om_equi.h'
*----------------------------------------------------------------------
C	EQUIV2
C		computes tridiagonal matrix

C	Input	Desc
C	NPTS	Number of points in vertical
C	Z	NPTS array of depths
C	NSQ	NPTS array of squared buoyancy frequency
C	GPINV	reduced gravity for bottom layer

C	Output	Desc
C	IPRIME	number of rows (and columns) in matrix
C	DEQ	Diagonal of tridiagonal matrix
C	UEQ	upper diagonal
C	L	lower diagonal

C	Internal arrays
C	LOC		mapping from reduced depth array to 
C			original depth array

      INTEGER LOC(MAXZVAL) 
      
C     ******************************************************************
C     *  computes tridiagonal matrix                                   *
C     *  free surface BC                                               *
C     ******************************************************************
C	DEQ has the negative of diagonal of the tridiagonal matrix [T]
C	where 
C		[T] G = - lambda G
C
C	and lambda is the eigenvalue.
      
C	UEQ is the upper diagonal
C	L   is the lower diagonal
      
      IF(GPINV.EQ.0D0)THEN
         NPTDM = NPTS
      ELSE
         NPTDM = NPTS - 1
      END IF
      
      DEQ(1) = GRAINV / DZ(2)
      UEQ(2) = DEQ(1)
C
      L(1)   = 0D0
      UEQ(1) = 0D0
      LOC(1) = 1
      IMIN   = 2
      IPRIME = IMIN
      IMAX = NPTDM - 1


  490 CONTINUE

      DO 500 I = IMIN, IMAX
C           cks for N = 0 (i.e. checks for zero buoyancy frequency)
         NCUR = NSQ(I)
         IF ( NCUR .LT. 1E-20 ) GO TO 510
C           calculates diagonal
         DEQ(IPRIME)= 2D0 / ( NCUR * DZ(I) * DZ(I+1) )
C           calculates off-diagonal
         L(IPRIME) = 2D0 / ( NCUR * DZ(I) * (DZ(I) + DZ(I+1)) )
         UEQ(IPRIME+1) = 2D0 / ( NCUR * DZ(I+1) * (DZ(I) + DZ(I+1)) )
         LOC(IPRIME) = I
         IPRIME = IPRIME + 1
  500 CONTINUE

      GOTO 540

C     *  NCUR = 0 -- eliminates pts from matrix                        *
  510 IMIN = I
C     *  searchs for nonzero NSQ                                       *
      DO 520 J = IMIN, IMAX
         IF ( NSQ(J) .GT. 1E-20 ) GO TO 530
  520 CONTINUE

C     ******************************************************************
C     *  mixed layer goes to bottom of profile                         *
C     ******************************************************************
      
      J = NPTDM
C     ******************************************************************
C     *  has found nonzero NSQ                                         *
C     *  the mode is linear over the range where NSQ is zero           *
C     *  ( from I-1 to J ).                                            *
C     *  Transforms (I-1)th equation                                   *
C     ******************************************************************
  530 DELTAZ = Z(J) - Z(I-1)
      DEQ(IPRIME-1) = DEQ(IPRIME-1) - UEQ(IPRIME) * (Z(J)-Z(I)) / DELTAZ
      IF ( J .EQ. NPTDM ) GOTO 540
      
      UEQ(IPRIME) = UEQ(IPRIME) * (Z(I)-Z(I-1)) / DELTAZ
C     *  calculates Jth equation and transforms                        *
      LOC(IPRIME) = J
      NCUR = NSQ(J)
      DEQ(IPRIME)   = 2.0 / ( NCUR * DZ(J) * DZ(J+1) )
      L(IPRIME)   = 2.0 / ( NCUR * DZ(J) * (DZ(J)+DZ(J+1)) )
      UEQ(IPRIME+1) = 2.0 / ( NCUR * DZ(J+1) * (DZ(J)+DZ(J+1)) )
C
      DEQ(IPRIME) = DEQ(IPRIME) - L(IPRIME) * (Z(J-1)-Z(I-1)) / DELTAZ
      L(IPRIME) =   L(IPRIME) * (Z(J)-Z(J-1)) / DELTAZ
C     ******************************************************************
C     *  returns to tridiagonal matrix calculating loop                *
C     ******************************************************************
      IPRIME = IPRIME + 1
      IMIN = J + 1
      IF ( IMIN .LE. IMAX ) GO TO 490
      
      
C     *  imposes bottom BC
C	* homogeneous if gpinv = 0.0
C	* otherwise there is an equivalent lower layer
C	* the equivalent lower layer has thickness dz(npts+1)
C	* and a layer buoyancy 1.0/gpinv
  540 CONTINUE
      L(IPRIME) = .5 * GPINV / DZ(NPTDM)
      DEQ(IPRIME) = .5 * GPINV * ( 1.0/DZ(NPTDM)
     *     + 1.0/DZ(NPTDM+1)   )
      LOC(IPRIME) = NPTDM

c      print*,dz(1:34)
c      print*,L(1:34)
c      print*,UEQ(1:34)
c      print*,DEQ(1:34)
c      stop
      
      return
      END
      
C	EQUIV3
C		computes eigenvalues

C	Inputs	Desc
C	IPRIME	number of rows (and columns) in matrix
C	DEQ	Diagonal of tridiagonal matrix
C	UEQ	upper diagonal
C	L	lower diagonal
C	NMODES	number of eigenvalues to be calculated

C	Output
C	EV	eigenvalues

*----------------------------------------------------------------------
      subroutine EQUIV3

      INCLUDE 'om_equi.h'
*----------------------------------------------------------------------
C	combines off-diagonals for eigenvalue routine

      IF(GPINV.EQ.0D0)THEN
         I1 = IPRIME - 1
      ELSE
         I1 = IPRIME
      END IF
      DO 550 I = 2, I1
         UEQ(I) = L(I) * UEQ(I)
  550 CONTINUE
      UEQ(1) = 0
C     sets pos. def. flag to no
      ISW   = 0
      IER   = 0
      
      
C  The DEQ array is equivalenced to the EV array.
C  Input to EQRT1S is from DEQ and UEQ.
C  Eigenvalue output goes to EV 
      
      CALL EQRT1S( DEQ, UEQ, I1, NMODES, ISW, IER )
      IF(IER.NE.0)THEN
         WRITE(0,*) 'IER=',IER,' in EQUIV3'
         STOP
      ENDIF
      
  570 CONTINUE
      
      return
      END
      
*----------------------------------------------------------------------
      subroutine EQUIV4(ZIS,ZTOP,VGRID,NZ,hstarp,
     *  CPH,LIGHT,SLVL,HVM,VVM,HCM,hvmsf,ahmix0,ahmixp)
      INCLUDE 'om_equi.h'
*----------------------------------------------------------------------
C	equiv4		- subroutine to calculate eigenvectors
C
C	The DO loop at line 700 does not compile properly
C	with optimization levels 2 or 3
C	(unless you include a write statement in the loop)
C
C	Input -  	Desc
C	NMODES	common	number of modes (includes barotropic)
C	EV	common	eigenvalues
C	NSQ	common	buoyancy frequency squared in rad/sec squared
C	TEMP1	common  temperature	
C	DZ	common	array of depth differences
C	GPINV	common  inverse reduced gravity for bottom layer
C
C	HVM	selected horizontal velocity mode values
C	VVM	selected vertical velocity mode values
C	HCM	selected heat content mode values
C
C	AHMIX0	meqsf	depth of T=20 degrees
C	AHMIXP	meqsf	vertical velocity mode value at z=ahmix0

      CHARACTER ZIS*(*)
      REAL CPH(*),LIGHT(*),SLVL(*), VGRID(NZ)
      REAL HVM(NZ,*),VVM(NZ,*),HCM(NZ,*)
      REAL ahmix0(1),ahmixp(*),hvmsf(*)
      
C	EMHOR	horizontal velocity mode
C	VR	vertical velocity mode (not in common)
C	HC	heat content mode	
      REAL VR(MAXZVAL),EMHOR(MAXZVAL),HC(MAXZVAL)
      real THTEMP(1)
      
      REAL*8 VRDP(MAXZVAL) , WGHT , D1SQ,D2SQ,DZ2,SUMA
C	variable to pass grids to LININT
      REAL VIN(MAXZVAL)
      
C	decides on grid for interpolated variables
      IF(ZIS(1:1).EQ.'T')THEN
C	temperature grid
         DO 100 IZ = 1 , NPTS
            VIN(IZ) = TEMP1(IZ)
  100    CONTINUE
         
      ELSE
C	depth grid
         DO 105 IZ = 1 , NPTS
            VIN(IZ) = Z(IZ)
  105    CONTINUE
      ENDIF
C	output grid is based on ZTOP DZR NZ
      ZVAL = ZTOP
      
      
C     ******************************************************************
C     *  Calculates eigenvectors and prints out                        *
C     *   loops over modes                                             *
C     ******************************************************************
      
C CALCULATION OF HSTARP IN THE LOOP - 9.1.86 CHRIS
      HSTARP=0.
      
      DO 1400 JMODE = 2, NMODES
         IMODE = JMODE - 1
C	JMODE points to all the modes,
C		i.e., it keeps track of the data in the arrays.
C	IMODE points to the baraclinic modes, only.
         
         PHASEV = SQRT( ABS( 1.0 / EV(JMODE) ) )
         CPH(IMODE) = PHASEV
C
C  NEW CODE, 04 Feb 1986
         
         IF(GPINV.NE.0D0)THEN
            VRDP(NPTS-2) = 1D0
            WGHT = GPINV*(1D0/DZ(NPTS-1) + 1D0/DZ(NPTS))-2D0*EV(JMODE) 
            VRDP(NPTS-1) = GPINV / (WGHT * DZ(NPTS-1) )
            VRDP(NPTS) = 0.0D0
            IMAX = NPTS - 2
         ELSE
            VRDP(NPTS-1) = 1D0
            VRDP(NPTS) = 0D0
            IMAX = NPTS - 1
         ENDIF
         IOLD = IMAX + 1
         DO 700 I = IMAX, 2, -1
            VRDP(I-1)  = VRDP(I) + ( DZ(I)/DZ(IOLD) )*(VRDP(I)-VRDP(IOLD) ) -
     *           EV(JMODE) * NSQ(I) * DZ(I) * ( (DZ(I)+DZ(IOLD))/2.0 ) *
     *           VRDP(I)
            IOLD = I
C	WRITE(0,*)I,VRDP(I-1),NSQ(I),DZ(I),DZ(I+1)
  700    CONTINUE
         
C     Normalizes
C
C	Calculate horizontal mode
         EMHOR(1) = -(VRDP(2)-VRDP(1)) / DZ(2)
         SUMA = EMHOR(1)**2 *  DZ(2)
         
         DO 990 I = 2, IMAX
  950       D1SQ = DZ(I) * DZ(I)
            D2SQ = DZ(I+1) * DZ(I+1)
            DZ2  = DZ(I) + DZ(I+1)
            EMHOR(I) = -( D1SQ*VRDP(I+1)-D2SQ*VRDP(I-1)+(D2SQ-D1SQ)*VRDP(I))
     W           / ( DZ(I) * DZ(I+1) * DZ2 )
            SUMA = SUMA + EMHOR(I)**2 * DZ2
  990    CONTINUE
         IF(GPINV.NE.0.0)THEN
            I = NPTS - 1
            EMHOR(I) = VRDP(I)*(1.0/DZ(NPTS) - EV(JMODE)/GPINV)
            SUMA = SUMA + EMHOR(I)**2 * DZ(I)
            I = NPTS
            EMHOR(I) = VRDP(I-1) / DZ(I)
            SUMA = SUMA + EMHOR(I)**2 * DZ(I)*2.
         ELSE
            I = NPTS
            EMHOR(I) = VRDP(I-1) / DZ(I)
            SUMA = SUMA + EMHOR(I)**2 * DZ(I)
         ENDIF
         SUMA = SUMA * 0.5
         FAC=SQRT(Z(NPTS)/SUMA)
         IF(EMHOR(1).LT.0.0)FAC=-FAC
         
C	Normalization
C	If we write the vertical velocity eigenfunction as G(z),
C	then the horizontal velocity eigenfunction is Gz(z).
C	The normalization is such that
C	integral from -D to 0 of [Gz(z)] squared = D
C	In a constant buoyancy frequency ocean that choice means 
C	the maximum amplitude is the square root of 2.
C	The normalization also means that G(z) has units of meters.
         
         DO 995 J=1, NPTS
            EMHOR(J)=FAC*EMHOR(J)
  995    CONTINUE
         XDISP = PHASEV/NBAR
C	normalizes vertical velocity/displacement mode as well
         DO 996 J = 1 , NPTS
            VR(J) = FAC * VRDP(J)
  996    CONTINUE
         
CC  CALCULATE LIGHTHILL COEFFICIENTS USING TRAPEZOIDAL RULE
CC  SUM=INTEGRAL OF EMHOR FROM Z=0 TO Z=-DML
C      SUM = 0.0
C
C      DO 1310 J = 1, NPTS
C         IF ( Z(J+1) .GT. DML ) GO TO 1320
C         SUM = SUM + DZ(J+1) * (EMHOR(J+1)+EMHOR(J)) / 2.0
C1310  CONTINUE
C      GO TO 1330
C
C1320  SUM = SUM + EMHOR(J) * ( DML-Z(J) )
C      SUM = SUM + ( (EMHOR(J+1)-EMHOR(J)) * (DML-Z(J))**2 ) /
C     * ( 2. * DZ(J+1) )
C
C1330  SUMCL = SUM / DML
C	CL is the coupling coef. (M):
C	DU/DT + ... = TAUX/CL
C      CL = Z(NPTS) / ( SUMCL * EMHOR(1) )
         
C	computes CL using just the surface value
         CL = Z(NPTS) / EMHOR(1)**2
         LIGHT(IMODE)=CL
c         print*,'lighthill coef (m)',imode,light(imode)
         HSTARP=HSTARP+1./CL
C	calculates SEALVL
         SLVL(IMODE) = EMHOR(1)/(EV(JMODE)*GRAV)
         
         HVMSF(IMODE) = EMHOR(1)
         
C	saves value at z0 = 50m (basic state profile depth)
CC	matches Steve's and Richard's TDEEP H input
         THTEMP(1) = 20.
c  really only need next line for one mode, oh well...
         CALL LININT (NPTS, TEMP1, Z, 1, THTEMP, AHMIX0, NLOW, NHIGH) 
c         print*,'calling LININT',AHMIX0
         CALL LININT (NPTS, TEMP1, VR, 1, THTEMP, AHMIXP(IMODE)
     *        ,NLOW, NHIGH)
c         print*,'calling LININT',VRDP(1),VRDP(2),VRDP(3),VRDP(4)
c         stop
         
         CALL LININT (NPTS , VIN , EMHOR ,NZ,VGRID, HVM(1,IMODE)
     *        ,NLOW , NHIGH )
         CALL LININT (NPTS , VIN , VR ,NZ,VGRID, VVM(1,IMODE)
     *        ,NLOW , NHIGH )

         
C	computes heat content
         SUM = 0.0
         HC(1) = SUM
         X = EMHOR(1)*TEMP1(1) *.5
         DO 1340 J = 2 , NPTS
            XL = X
            X = EMHOR(J) * TEMP1(J) * .5
            SUM = SUM + (X+XL)*DZ(J)
            HC(J) = SUM
 1340    CONTINUE
         CALL LININT(NPTS,VIN,HC, NZ,VGRID,HCM(1,IMODE), NLOW,NHIGH)
         
 1400 CONTINUE
      TLIGHT = 1./HSTARP
C	WRITE(0,*) DML,TLIGHT
c      print*,'TLIGHT = 1/HSTARP',TLIGHT

      return
      END
      
C	linear interpolation to new grid
C	inputs
C	NT	number of points in the input profile
C	X	array of x
C	Y	array of y
C	NOUT
C	XOUT
C	outputs
C	NLOW	index of lowest XOUT within X
C	NHIGH	index of highest XOUT within X
C	YOUT	Y regridded to XOUT
C	
	subroutine LININT(NT,X,Y,NOUT,XOUT,YOUT,NLOW,NHIGH)
	REAL X(NT) , Y(NT) , XOUT(NOUT) , YOUT(NOUT)
C	initializes 
	XBINL = X(1)
	XBINH = X(NT)
	DELTAX = XBINH-XBINL
	DXOUT = XOUT(NOUT)-XOUT(1)
		IF((DXOUT*DELTAX).GT.0.0)THEN
C	both ascending or descending in X and XOUT
	IMIN = 1
	IMAX = NOUT
	ISTP = 1
				ELSE
C	opposite
	IMIN = NOUT
	IMAX = 1
	ISTP = -1
					ENDIF
	IFLIP = - ISTP

C	two versions -- one ascending in x, one descending
					IF(DELTAX.GT.0.0)THEN

	DO 20 I = IMIN , IMAX , ISTP
	IF(XBINL.LE.XOUT(I))GOTO 25
	YOUT(I) = Y(1)
20	CONTINUE
	I = 0
25	NLOW = I
C	finds NHIGH
	DO 30 I = IMAX , IMIN , IFLIP
	IF(XBINH.GE.XOUT(I))GOTO 35
	YOUT(I) = Y(NT)
30	CONTINUE
	I = 0
35	NHIGH = I
C	interpolates
	XBINH = X(2)
	IPTH = 2
	IPTL = 1
	DELTAX = XBINH-XBINL
	IF(NLOW.EQ.0)GOTO 8888
	IF(NHIGH.EQ.0)GOTO 8888
      DO 100 ICNT = NLOW , NHIGH , ISTP
	XVAL = XOUT(ICNT)
50	CONTINUE
			      IF(XVAL.GT.XBINH)THEN
	IPTL = IPTH
	IPTH = IPTH + 1
	XBINL = XBINH
	XBINH = X(IPTH)
	DELTAX = XBINH-XBINL
	GOTO 50
					     END IF

      PERCENT = ( XVAL - XBINL )  / DELTAX
      YOUT(ICNT) = Y(IPTL) + PERCENT * ( Y(IPTH) - Y(IPTL) )

  100 CONTINUE
							ELSE
C	descending
	DO 120 I = IMIN , IMAX , ISTP
	IF(XBINL.GE.XOUT(I))GOTO 125
	YOUT(I) = Y(1)
120	CONTINUE
	I = 0
125	NLOW = I
C	finds NHIGH
	DO 130 I = IMAX , IMIN , IFLIP
	IF(XBINH.LE.XOUT(I))GOTO 135
	YOUT(I) = Y(NT)
130	CONTINUE
	I = 0
135	NHIGH = I
C	interpolates
	XBINH = X(2)
	IPTH = 2
	IPTL = 1
	DELTAX = XBINH-XBINL
	IF(NLOW.EQ.0)GOTO 8888
	IF(NHIGH.EQ.0)GOTO 8888
      DO 200 ICNT = NLOW , NHIGH , ISTP
	XVAL = XOUT(ICNT)
150	CONTINUE
			      IF(XVAL.LT.XBINH)THEN
	IPTL = IPTH
	IPTH = IPTH + 1
	XBINL = XBINH
	XBINH = X(IPTH)
	DELTAX = XBINH-XBINL
	GOTO 150
					     END IF

      PERCENT = ( XVAL - XBINL )  / DELTAX
      YOUT(ICNT) = Y(IPTL) + PERCENT * ( Y(IPTH) - Y(IPTL) )
200   CONTINUE
							END IF

8888	CONTINUE
C	patch so that NLOW is always the lower index if the 
C		range is not empty
	IF(ISTP.LT.0)THEN
	I = NLOW
	NLOW = NHIGH
	NHIGH = I
		ENDIF
      return
      END

!$Source: /home/benno/Models/LinearEquator/RCS/eqrt1s.f,v $
!$Author: benno $
!$Date: 1995/03/09 12:12:38 $
!$Revision: 1.1 $
!$State: Exp $
      subroutine EQRT1S (D,E2,N,M,ISW,IER)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,ISW,IER
      REAL*8               D(1001),E2(1001)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II,I,JJ,J,K1,K
      REAL*8               DELTA,DLAM,EP,ERR,F,P,QP,Q,RLDLP,R,S,TOT
C                                  RLDLP = MACHINE PRECISION
      DATA               RLDLP/1.192093E-07/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DLAM = 0.0E0
      ERR = 0.0E0
      S = 0.0E0
C                                  LOOK FOR SMALL SUB-DIAGONAL ENTRIES
C                                  DEFINE INITIAL SHIFT FROM LOWER
C                                  GERSCHGORIN BOUND.
      TOT = D(1)
      Q = 0.0E0
      J = 0
      DO 15 I=1,N
         P = Q
         IF (I.EQ.1) GO TO 5
         IF (P.GT.RLDLP*(ABS(D(I))+ABS(D(I-1)))) GO TO 10
    5    E2(I) = 0.0E0
C                                  COUNT IF E2(I) HAS UNDERFLOWED
   10    IF (E2(I).EQ.0.0E0) J = J+1
         Q = 0.0E0
         IF (I.NE.N) Q = SQRT(E2(I+1))
         TOT = DMIN1(D(I)-P-Q,TOT)
   15 CONTINUE
      IF (ISW.EQ.1.AND.TOT.LT.0.0E0) GO TO 25
      DO 20 I=1,N
   20 D(I) = D(I)-TOT
      GO TO 30
   25 TOT = 0.0E0
   30 DO 90 K=1,M
C                                  NEXT QR TRANSFORMATION
   35    TOT = TOT+S
         DELTA = D(N)-S
         I = N
         F = ABS(RLDLP*TOT)
         IF (DLAM.LT.F) DLAM = F
         IF (DELTA.GT.DLAM) GO TO 40
         IF (DELTA.GE.(-DLAM)) GO TO 75
         IER = 130
         GO TO 9000
C                                  REPLACE SMALL SUB-DIAGONAL SQUARES
C                                  BY ZERO TO REDUCE THE INCIDENCE OF
C                                  UNDERFLOWS
   40    IF (K.EQ.N) GO TO 50
         K1 = K+1
         DO 45 J=K1,N
            IF (E2(J).LE.(RLDLP*(D(J)+D(J-1)))**2) E2(J) = 0.0E0
   45    CONTINUE
   50    F = E2(N)/DELTA
         QP = DELTA+F
         P = 1.0E0
         IF (K.EQ.N) GO TO 65
         K1 = N-K
         DO 60 II=1,K1
            I = N-II
            Q = D(I)-S-F
            R = Q/QP
            P = P*R+1.0E0
            EP = F*R
            D(I+1) = QP+EP
            DELTA = Q-EP
            IF (DELTA.GT.DLAM) GO TO 55
            IF (DELTA.GE.(-DLAM)) GO TO 75
            IER = 130
            GO TO 9000
   55       F = E2(I)/Q
            QP = DELTA+F
            E2(I+1) = QP*EP
   60    CONTINUE
   65    D(K) = QP
         S = QP/P
         IF (TOT+S.GT.TOT) GO TO 35
         IER = 33
         E2(1) = K
C                                  SET ERROR -- IRREGULAR END
C                                  DEFLATE MINIMUM DIAGONAL ELEMENT
         S = 0.0E0
         DELTA = QP
         DO 70 J=K,N
            IF (D(J).GT.DELTA) GO TO 70
            I = J
            DELTA = D(J)
   70    CONTINUE
C                                  CONVERGENCE
   75    IF (I.LT.N) E2(I+1) = E2(I)*F/QP
         IF (I.EQ.K) GO TO 85
         K1 = I-K
         DO 80 JJ=1,K1
            J = I-JJ
            D(J+1) = D(J)-S
            E2(J+1) = E2(J)
   80    CONTINUE
   85    D(K) = TOT
         ERR = ERR+ABS(DELTA)
         E2(K) = ERR
   90 CONTINUE
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE
      CALL UERTST (IER,'EQRT1S')
 9005 return
      END
C   IMSL ROUTINE NAME   - UERTST
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - ELXSI/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
C
C   USAGE               - CALL UERTST (IER,NAME)
C
C   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
C                           IER = I+J WHERE
C                             I = 128 IMPLIES TERMINAL ERROR MESSAGE,
C                             I =  64 IMPLIES WARNING WITH FIX MESSAGE,
C                             I =  32 IMPLIES WARNING MESSAGE.
C                             J = ERROR CODE RELEVANT TO CALLING
C                                 ROUTINE.
C                NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING
C                           THE NAME OF THE CALLING ROUTINE. (INPUT)
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - UGETIO,USPKD
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
C                TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
C                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
C                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
C                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
C                UGETIO AS FOLLOWS..
C                                NIN = 0
C                                NOUT = NEW OUTPUT UNIT NUMBER
C                                CALL UGETIO(3,NIN,NOUT)
C                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      subroutine UERTST (IER,NAME)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      CHARACTER          NAME*(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEQDF,IOUNIT,LEVEL,LEVOLD,NIN,NMTB
      CHARACTER          IEQ,NAMEQ(6),NAMSET(6),NAMUPK(6)
      DATA               NAMSET/'U','E','R','S','E','T'/
      DATA               NAMEQ/6*' '/
      DATA               LEVEL/4/,IEQDF/0/,IEQ/'='/
C                                  UNPACK NAME INTO NAMUPK
C                                  FIRST EXECUTABLE STATEMENT
      CALL USPKD (NAME,6,NAMUPK,NMTB)
C                                  GET OUTPUT UNIT NUMBER
      CALL UGETIO(1,NIN,IOUNIT)
C                                  CHECK IER
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
C                                  PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
C                                  PRINT WARNING WITH FIX MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
C                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
      GO TO 30
   15 CONTINUE
C                                  CHECK FOR UERSET CALL
      DO 20 I=1,6
         IF (NAMUPK(I).NE.NAMSET(I)) GO TO 25
   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
C                                  PRINT NON-DEFINED MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK
   30 IEQDF = 0
      return
   35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   40 FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
   50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     1       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
C
C                                  SAVE P FOR P = R CASE
C                                    P IS THE PAGE NAMUPK
C                                    R IS THE ROUTINE NAMUPK
   55 IEQDF = 1
      DO 60 I=1,6
   60 NAMEQ(I) = NAMUPK(I)
   65 return
      END
C   IMSL ROUTINE NAME   - UGETIO
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - ELXSI/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1981
C
C   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
C                           VALUES FOR INPUT AND OUTPUT UNIT
C                           IDENTIFIERS.
C
C   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
C
C   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
C                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
C                           UNIT IDENTIFIER VALUES ARE returnED IN NIN
C                           AND NOUT, RESPECTIVELY.
C                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
C                           RESET FOR SUBSEQUENT USE.
C                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
C                           RESET FOR SUBSEQUENT USE.
C                NIN    - INPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
C                NOUT   - OUTPUT UNIT IDENTIFIER.
C                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
C                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
C                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
C                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
C                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      subroutine UGETIO(IOPT,NIN,NOUT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/5/,NOUTD/6/
C                                  FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
   10 NOUTD = NOUT
 9005 return
      END
C   IMSL ROUTINE NAME   - USPKD
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - ELXSI/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE
C                           CHARACTER STRING ARGUMENTS
C
C   USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C
C   ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)
C                NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.
C                UNPAKD - CHARACTER ARRAY TO RECEIVE THE UNPACKED
C                         REPRESENTATION OF THE STRING. (OUTPUT)
C                NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)
C
C   PRECISION/HARDWARE  - SINGLE/ALL
C
C   REQD. IMSL ROUTINES - NONE
C
C   REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO A CHARACTER ARRAY
C                IN (A1) FORMAT.
C            2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF
C                THAT ARE IGNORED.
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
      subroutine USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,NCHARS,NCHMTB
C
      CHARACTER          UNPAKD(1),IBLANK
      CHARACTER*(1)      PACKED(1)
      DATA               IBLANK /' '/
C                                  INITIALIZE NCHMTB
      NCHMTB = 0
C                                  return IF NCHARS IS LE ZERO
      IF(NCHARS.LE.0) return
C                                  SET NC=NUMBER OF CHARS TO BE DECODED
      NC = MIN0 (129,NCHARS)
      DO 5 I=1,NC
         UNPAKD(I) = PACKED(I)
    5 CONTINUE
  150 FORMAT (129A1)
C                                  CHECK UNPAKD ARRAY AND SET NCHMTB
C                                  BASED ON TRAILING BLANKS FOUND
      DO 200 N = 1,NC
         NN = NC - N + 1
         IF(UNPAKD(NN) .NE. IBLANK) GO TO 210
  200 CONTINUE
  210 NCHMTB = NN
      return
      END
