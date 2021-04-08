!----------------------------------------------------------------------
      subroutine VMODESSUB(NIN,PIN,ANBV2,TEMP,NUMMOD,ZCUT,ZIS,ZTOP,VGRID,
     *     NZ,hstarp,CPH,LIGHT,SLVL,HVM,VVM,HCM,HVMSF,ahmix0,ahmixp)
      
      INCLUDE 'om_equi.h'
!----------------------------------------------------------------------
! NAME
!	vmodessub(9)		- subroutine to calculate vertical modes
! SYNOPSIS
!	subroutine VMODESSUB(NIN,PIN,ANBV2,TEMP,NUMMOD,ZCUT,
!     *  ZIS,ZTOP,VGRID,NZ,hstarp,CPH,LIGHT,SLVL,HVM,VVM,HCM)
! Inputs
!	NIN	[]	number of points in input profile
!	PIN	[P>	input pressures (or depths (m))
!	ANBV2	[ANBV2>	input BV frequency squared (rad/s)squared
!	TEMP	[TEMP>	input temperature profile
!	NUMMOD	()	number of modes to calculate
!	ZCUT		equivalent bottom layer depth parameter
!			has no effect if greater than water col depth
!	ZIS	char	vertical coordinate name
!			either Depth or Temperature
!	ZTOP		smallest coordinate to interpolate output to
!	VGRID		output grid
!	NZ	{}	number of points in output grid
! Outputs
!	hstarp   	1./(total Lighthill coefficient over NUMMOD modes)
!	CPH	(CPH>	phase speeds (m/s)
!	LIGHT	(LIGH>	Lighthill coefficents (1/m)
!	SLVL	(SLVL>	Sea Level (m)
!	HVM	{HVM)	Horizontal Velocity Modes
!	VVM	{VVM)	Vertical Velocity Modes
!	HCM	{HCM)	Heat Content Modes
!     HVMSF           surface value of HVM
!	ahmix0		depth of 20degree isotherm
!	ahmixp		value of vertical velocity mode at 20degree isotherm
!
!	all arguments are real*4
!	double precision imsl library required
!
! DESCRIPTION
!	This program uses a tridiagonal matrix solver to calculate dynamical
!	modes from a buoyancy frequency profile.  The problem is solved with
!	a free surface boundary condition.
! BUGS
!	Clearly this is an awkward subroutine, with way too many arguments.
! AUTHOR
!       Benno provided this version of the equiv* subroutines
!----------------------------------------------------------------------

      REAL	PIN(NIN),ANBV2(NIN),TEMP(NIN)
      REAL	CPH(NUMMOD),LIGHT(NUMMOD),SLVL(NUMMOD)
      REAL	HVM(NZ,NUMMOD),VVM(NZ,NUMMOD),HCM(NZ,NUMMOD)
      REAL	ahmix0(NUMMOD),ahmixp(NUMMOD),hvmsf(nummod)
      
      REAL	ZCUT,ZTOP,VGRID(NZ),hstarp
      
      CHARACTER*(*) ZIS
      CHARACTER*256	FNAME
      CHARACTER*80	LINE
      
!	altered 27 April 1987 to include from local directory
!	19 May 1987 various modifications including equivalent bottom layer,
!		double precision eigenvalue/vector calculation, temperature
!		coordinates, and new printout
!	altered 15 AUG 1997 - Naomi- arbitrary output grid, vgrid

      GRAV   = 9.78049
      GRAINV = 1.0 / GRAV
      
! DML is the mixed layer depth. disabled!!!!!!!!!!!!
!      DML = HSFC
! 	NOUT is the COMMON block copy of NZ
      NOUT = NZ
! NMODES includes the barotropic mode (#1) which will be discarded at the end.
! NUMMOD includes only baroclinic modes, (#2 to #n).
      NMODES = NUMMOD + 1
      
!  Read the water profile and...
      CALL EQUIV1(ZCUT,ZIS,ZTOP,VGRID,NZ,NIN,PIN,ANBV2,TEMP)
      
!  Compute the tridiagonal matrix
      CALL EQUIV2
      
!  Call the eigenvalue routine
      CALL EQUIV3
      
!  calculate eigenvectors, fns of eigenvalues
      CALL EQUIV4(ZIS,ZTOP,VGRID,NZ,hstarp,
     *     CPH,LIGHT,SLVL,HVM,VVM,HCM,hvmsf,ahmix0,ahmixp)
!	HVM	{HVM)	selected horizontal velocity mode values
!	VVM	{VVM)	selected vertical velocity mode values
!	HCM	{HCM)	selected heat content mode values
      
      return
      END
      
!----------------------------------------------------------------------
      subroutine MODEPA(nx,ny,nmodes,npt,time_mul,delt,delx,dely,xm,ym,
     *            X_MIN,X_MAX,Y_MIN,Y_MAX,cph,phasev,cl,sealvl,dtm,A_rayl,rayl,
     *            A_lapl,rlap,hxm,hym,fcor,fc,hscl,uscl,dscl,ahmixp,
     *            scale_L)
      implicit none
!----------------------------------------------------------------------
! manpage modepa
! NAME
!	modepa(9)[vmodessub]	- subroutine to calculate scales for modes
! SYNOPSIS
!	subroutine MODEPA([grid_stuff],PHASEV,SEALVL,CL,HSCL,USCL,DSCL)
! Inputs
!       nx              number of zonal grid points 
!       ny              number of meridional grid points 
!       nmodes          number of vertical modes
!       time_mul        1: reset delt = time_mul*delx 
!       delt            time step (s)
!       A_rayl          Rayleigh friction coefficient
!       A_lapl          Laplacian diffusion coefficient
!       X_MIN,X_MAX     zonal limits (degrees)
!       Y_MIN,Y_MAX     meridional limits (degrees)
!	phasev		phase velocity (m/s)
!	sealvl		sea level coefficient
!	cl 		Lighthill coupling coef.:	du/dt + ... = taux/CL
!	ahmixp 		vertical velocity eigenfunction at z(temp=20)
! Outputs
!       delt            modified time step (s)
!       delx            zonal grid spacing (m)
!       dely            merid grid spacing (m)
!       xm              zonal grid (m)
!       ym              merid grid (m)
!       dtm             scaled time
!       rayl            scaled rayleigh friction
!       rlap            scaled laplacian coefficient
!       hxm             scaled zonal grid spacing
!       hym             scaled merid grid spacing
!       f!              scaled coriolis parameter
!       fcor            unscaled coriolis parameter (1/s)
!	hscl		scaling for pressure (m)
!	uscl		scaling for velocity (m/s)
!	dscl		scaling for vertical displacement (m)
!	ahmixp 		scaled ahmixp for tdeep
! DESCRIPTION
!	modepa is a subroutine to calculate mode scales from the information
!	that vmodessub returns.
! BUGS
!	It has not been tested yet
! AUTHOR
!       Benno provided the skeleton of this subroutine
!----------------------------------------------------------------------

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

!               Rayleigh friction, scaled for each mode
         rayl(mode) = rt * A_rayl * (phasev(1)/phasev(mode))**2
         rlap(mode) = rt * A_lapl / rl**2

!       	scale for displacement (isopycnal height) (m)
         do i = 1, nx*ny
            dscl(mode,i)= grainv* rl/ (cl(mode,i)* sealvl(mode,i) )
            
! precompute coefs for tdeep:
            ahmixp(mode,i) =-dscl(mode,i) * ahmixp(mode,i) 
         enddo
            
!	        scale for velocity  (u,v) (m/s)
         uscl(mode) = phasev(mode)
      
!	        scale for pressure perturbations (h) (m)
         hscl(mode) = grainv * (phasev(mode)**2)

!	        scaled time step
         dtm(mode) = delt/rt

!	        scaled zonal grid spacing
         hxm(mode) = delx/rl

!	        scaled meridional grid spacing
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
      
!----------------------------------------------------------------------
      subroutine single_mode(nx,ny,time_mul,delt,delx,dely,xm,ym,
     *          X_MIN,X_MAX,Y_MIN,Y_MAX,phasev,H,H1,hvm,hvmsf,dtm,A_rayl,rayl
     *          ,hxm,hym,fcor,fc,dscl,hscl,uscl,ahmix0,ahmixp,hstarp,area,
     *          scale_L)
      implicit none
!----------------------------------------------------------------------
! Inputs
!       nx              number of zonal grid points 
!       ny              number of meridional grid points 
!       time_mul        1: reset delt = time_mul*delx 
!       delt            time step (s)
!       A_rayl          Rayleigh friction coefficient
!       X_MIN,X_MAX     zonal limits (degrees)
!       Y_MIN,Y_MAX     meridional limits (degrees)
!	phasev		phase velocity (m/s)
!	H		depth of active ocean
!	H1		depth of surface layer
! Outputs
!       delt            modified time step (s)
!       delx            zonal grid spacing (m)
!       dely            merid grid spacing (m)
!       xm              zonal grid (m)
!       ym              merid grid (m)
!       dtm             scaled time
!       rayl            scaled rayleigh friction
!       hxm             scaled zonal grid spacing
!       hym             scaled merid grid spacing
!       f!              scaled coriolis parameter
!       fcor            unscaled coriolis parameter (1/s)
!	hscl		scaling for pressure (m)
!	uscl		scaling for velocity (m/s)
!	hvm		eigenfunction reset to 1.
!	hvmsf		surface value eigenfunction reset to 1.
!	ahmix0 		rescaled ahmix0
!	ahmixp 		rescaled ahmixp for tdeep
!----------------------------------------------------------------------

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
      
!	        scale for velocity  (u,v) (m/s)
      uscl = phasev
      
!	        scale for pressure perturbations (h) (m)
      hscl = H
      hvm(1) = 1.
      hvmsf(1)= 1.

      ahmix0(1) = H
      ahmixp(1) = hscl
      
!	        scaled time step
      dtm = delt/rt
      
!	        scaled zonal grid spacing
      hxm = delx/rl
      
!	        scaled meridional grid spacing
      hym = dely/rl

      area = rt/hscl/delx/dely
      
!         print*,'mode',mode,phasev(mode),dtm(mode),hxm(mode),hym(mode)

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
      
!----------------------------------------------------------------------

      subroutine EQUIV1(ZCUT,ZIS,ZTOP,VGRID,NZ,NPTSI,PIN,ANBV2,TEMP)
      
      INCLUDE 'om_equi.h'
!----------------------------------------------------------------------
!	equiv1[vmodes]	- profile processor for vmodes
! Inputs
!	NPTS	[]	number of points in input profile
!	PIN	[P>	input pressures
!	ANBV2	[ANBV2>	input BV frequency squared (rad/s)squared
!	TEMP	[TEMP>	input temperature profile
! Outputs
!	NBAR	common	Vertically averaged buoyancy frequency in rad/sec
!	Z	common  NPTS array of depths
!	NSQ	common  NPTS array of squared buoyancy frequency
      
      CHARACTER ZIS*(*)
      REAL PIN(NPTSI),ANBV2(NPTSI),TEMP(NPTSI), VGRID(NZ)
      
!	Internal arrays
!	EDENS1	equivalent density
!	HC0	basic state heat content
      
      REAL EDENS1(MAXZVAL),HC0(MAXZVAL)
      
      INTEGER	HEQLUN
!	variables to pass grids to LININT
      REAL VIN(MAXZVAL)
      
      NPTS = NPTSI

      DO 99 J = 1 , NPTS
         Z(J) = PIN(J)
         NSQ(J) = AMAX1(0.,ANBV2(J))
         TEMP1(J) = TEMP(J)
   99 CONTINUE
      
!     *  checks whether first point in profile is at Z=0   *
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
      
!	computes depth differences
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
      
!	computes equivalent density (integrated buoyancy frequency)
      SUM = 0.
      EDENS1(1) = 1.
      DO 500 I = 2,NPTS
         SUM = SUM + (NSQ(I)+NSQ(I-1))*DZ(I)*.5*GRAINV
         EDENS1(I) = EXP(SUM)
  500 CONTINUE
      
!	computes basic state heat content
      SUM = 0.
      HC0(1) = SUM
      DO 505 I = 2 , NPTS
         SUM = SUM + (TEMP1(I)+TEMP1(I-1))*DZ(I)*.5
         HC0(I) = SUM
  505 CONTINUE
      
!	uses full profile for TEMP and EDENS1 interpolation
!	even if profile gets cut
      
!	decides on grid for interpolated variables
      IF(ZIS(1:1).EQ.'T')THEN
!	temperature grid
         DO 510 IZ = 1 , NPTS
            VIN(IZ) = TEMP1(IZ)
  510    CONTINUE
      ELSE
!	depth grid
         DO 520 IZ = 1 , NPTS
            VIN(IZ) = Z(IZ)
  520    CONTINUE
      ENDIF
!	output grid is based on ZTOP and ZBOT
      ZVAL = ZTOP
      
!-----------------------------------------------------------------------
!	modifies profile with an equivalent bottom layer if requested
!	if ZCUT < Z(NPTS) then creates layer
      
      IF(ZCUT.GT.Z(NPTS))THEN
         GPINV = 0.0
         DZ(NPTS+1) = 1
      ELSE
!	rounds ZCUT to nearest gridpoint
         CALL BINSERCH ( ZCUT, Z, NPTS, IPT )
         IF(  (ZCUT-Z(IPT)) .GT. (Z(IPT+1)-ZCUT) )IPT = IPT + 1
         ZCUT = Z(IPT)
         
!	 computes inverse reduced gravity
         SUM = 0.0
         IMIN = IPT + 1
         IL = IPT
         DO 600 I = IMIN , NPTS
            SUM = SUM + (NSQ(I) + NSQ(IL))*DZ(I)
            IL = I
  600    CONTINUE
         GPINV = 2.0 / SUM
         
!	patchs bottom of profile
         DZ(IPT+1) 	= Z(NPTS) - ZCUT
         Z(IPT+1) 	= Z(NPTS)
         NSQ(IPT+1)	= NSQ(NPTS)
         TEMP1(IPT+1) 	= TEMP1(NPTS)
         
!	remaining patches unecessary since arrays are local to EQUIV1
         EDENS1(IPT+1)	= EDENS1(NPTS)
         HC0(IPT+1)	= HC0(NPTS)
         
         NPTS = IPT + 1
      END IF

      return
      END
      
!----------------------------------------------------------------------
      subroutine BINSERCH ( ZVAL, Z, NPTS, IPT )
!----------------------------------------------------------------------
!
!  03 June 1986
!  Search the Z array and get an IPT value.
!  We assume the Z array is sorted in ascending order.
!  If ZVAL is not between Z(1) and Z(NPTS) then IPT = 0.
!


      DIMENSION Z(NPTS)
      
!  In range?
      IF ( ZVAL .LT. Z(1) .OR. ZVAL .GT. Z(NPTS) ) THEN
         IPT = 0
         return
      END IF
      
! Start and end pointers, which will be 'squeezed together'.
      IS = 1
      IE = NPTS
      
! A 'WHILE' statement would be very nice at this point.
   40 CONTINUE
      
! The pointer and the spread/difference between the bounds
      IPT = ( IS + IE ) / 2
      IDIF = IE - IS
      
! A direct hit!
      IF ( ZVAL .EQ. Z(IPT) ) return
      
! Our bounds can't be squeezed any more.
      IF ( IDIF .EQ. 1 ) THEN
         IPT = IS
         return
      END IF
      
! Determine which bound to 'move in'.
      IF ( ZVAL .LT. Z(IPT) ) THEN
         IE = IPT
      ELSE IF ( ZVAL .GT. Z(IPT) ) THEN
         IS = IPT
      END IF
      
! Do it again.
      GO TO 40
      
      END

!----------------------------------------------------------------------
      subroutine EQUIV2
      INCLUDE 'om_equi.h'
!----------------------------------------------------------------------
!	EQUIV2
!		computes tridiagonal matrix
!
!	Input	Desc
!	NPTS	Number of points in vertical
!	Z	NPTS array of depths
!	NSQ	NPTS array of squared buoyancy frequency
!	GPINV	reduced gravity for bottom layer
!
!	Output	Desc
!	IPRIME	number of rows (and columns) in matrix
!	DEQ	Diagonal of tridiagonal matrix
!	UEQ	upper diagonal
!	L	lower diagonal
!
!	Internal arrays
!	LOC		mapping from reduced depth array to 
!			original depth array

      INTEGER LOC(MAXZVAL) 
      
!     ******************************************************************
!     *  computes tridiagonal matrix                                   *
!     *  free surface B!                                               *
!     ******************************************************************
!	DEQ has the negative of diagonal of the tridiagonal matrix [T]
!	where 
!		[T] G = - lambda G
!
!	and lambda is the eigenvalue.
!      
!	UEQ is the upper diagonal
!	L   is the lower diagonal
      
      IF(GPINV.EQ.0D0)THEN
         NPTDM = NPTS
      ELSE
         NPTDM = NPTS - 1
      END IF
      
      DEQ(1) = GRAINV / DZ(2)
      UEQ(2) = DEQ(1)
!
      L(1)   = 0D0
      UEQ(1) = 0D0
      LOC(1) = 1
      IMIN   = 2
      IPRIME = IMIN
      IMAX = NPTDM - 1


  490 CONTINUE

      DO 500 I = IMIN, IMAX
!           cks for N = 0 (i.e. checks for zero buoyancy frequency)
         NCUR = NSQ(I)
         IF ( NCUR .LT. 1E-20 ) GO TO 510
!           calculates diagonal
         DEQ(IPRIME)= 2D0 / ( NCUR * DZ(I) * DZ(I+1) )
!           calculates off-diagonal
         L(IPRIME) = 2D0 / ( NCUR * DZ(I) * (DZ(I) + DZ(I+1)) )
         UEQ(IPRIME+1) = 2D0 / ( NCUR * DZ(I+1) * (DZ(I) + DZ(I+1)) )
         LOC(IPRIME) = I
         IPRIME = IPRIME + 1
  500 CONTINUE

      GOTO 540

!     *  NCUR = 0 -- eliminates pts from matrix                        *
  510 IMIN = I
!     *  searchs for nonzero NSQ                                       *
      DO 520 J = IMIN, IMAX
         IF ( NSQ(J) .GT. 1E-20 ) GO TO 530
  520 CONTINUE

!     ******************************************************************
!     *  mixed layer goes to bottom of profile                         *
!     ******************************************************************
      
      J = NPTDM
!     ******************************************************************
!     *  has found nonzero NSQ                                         *
!     *  the mode is linear over the range where NSQ is zero           *
!     *  ( from I-1 to J ).                                            *
!     *  Transforms (I-1)th equation                                   *
!     ******************************************************************
  530 DELTAZ = Z(J) - Z(I-1)
      DEQ(IPRIME-1) = DEQ(IPRIME-1) - UEQ(IPRIME) * (Z(J)-Z(I)) / DELTAZ
      IF ( J .EQ. NPTDM ) GOTO 540
      
      UEQ(IPRIME) = UEQ(IPRIME) * (Z(I)-Z(I-1)) / DELTAZ
!     *  calculates Jth equation and transforms                        *
      LOC(IPRIME) = J
      NCUR = NSQ(J)
      DEQ(IPRIME)   = 2.0 / ( NCUR * DZ(J) * DZ(J+1) )
      L(IPRIME)   = 2.0 / ( NCUR * DZ(J) * (DZ(J)+DZ(J+1)) )
      UEQ(IPRIME+1) = 2.0 / ( NCUR * DZ(J+1) * (DZ(J)+DZ(J+1)) )
!
      DEQ(IPRIME) = DEQ(IPRIME) - L(IPRIME) * (Z(J-1)-Z(I-1)) / DELTAZ
      L(IPRIME) =   L(IPRIME) * (Z(J)-Z(J-1)) / DELTAZ
!     ******************************************************************
!     *  returns to tridiagonal matrix calculating loop                *
!     ******************************************************************
      IPRIME = IPRIME + 1
      IMIN = J + 1
      IF ( IMIN .LE. IMAX ) GO TO 490
      
      
!     *  imposes bottom BC
!	* homogeneous if gpinv = 0.0
!	* otherwise there is an equivalent lower layer
!	* the equivalent lower layer has thickness dz(npts+1)
!	* and a layer buoyancy 1.0/gpinv
  540 CONTINUE
      L(IPRIME) = .5 * GPINV / DZ(NPTDM)
      DEQ(IPRIME) = .5 * GPINV * ( 1.0/DZ(NPTDM)
     *     + 1.0/DZ(NPTDM+1)   )
      LOC(IPRIME) = NPTDM

!      print*,dz(1:34)
!      print*,L(1:34)
!      print*,UEQ(1:34)
!      print*,DEQ(1:34)
!      stop
      
      return
      END
      
!	EQUIV3
!		computes eigenvalues
!
!	Inputs	Desc
!	IPRIME	number of rows (and columns) in matrix
!	DEQ	Diagonal of tridiagonal matrix
!	UEQ	upper diagonal
!	L	lower diagonal
!	NMODES	number of eigenvalues to be calculated

!	Output
!	EV	eigenvalues

!----------------------------------------------------------------------
      subroutine EQUIV3

      INCLUDE 'om_equi.h'
!----------------------------------------------------------------------
!	combines off-diagonals for eigenvalue routine

      IF(GPINV.EQ.0D0)THEN
         I1 = IPRIME - 1
      ELSE
         I1 = IPRIME
      END IF
      DO 550 I = 2, I1
         UEQ(I) = L(I) * UEQ(I)
  550 CONTINUE
      UEQ(1) = 0
!     sets pos. def. flag to no
      ISW   = 0
      IER   = 0
      
      
!  The DEQ array is equivalenced to the EV array.
!  Input to EQRT1S is from DEQ and UEQ.
!  Eigenvalue output goes to EV 
      
      CALL EQRT1S( DEQ, UEQ, I1, NMODES, ISW, IER )
      IF(IER.NE.0)THEN
         WRITE(0,*) 'IER=',IER,' in EQUIV3'
         STOP
      ENDIF
      
  570 CONTINUE
      
      return
      END
      
!----------------------------------------------------------------------
      subroutine EQUIV4(ZIS,ZTOP,VGRID,NZ,hstarp,
     *  CPH,LIGHT,SLVL,HVM,VVM,HCM,hvmsf,ahmix0,ahmixp)
      INCLUDE 'om_equi.h'
!----------------------------------------------------------------------
!	equiv4		- subroutine to calculate eigenvectors
!
!	The DO loop at line 700 does not compile properly
!	with optimization levels 2 or 3
!	(unless you include a write statement in the loop)
!
!	Input -  	Desc
!	NMODES	common	number of modes (includes barotropic)
!	EV	common	eigenvalues
!	NSQ	common	buoyancy frequency squared in rad/sec squared
!	TEMP1	common  temperature	
!	DZ	common	array of depth differences
!	GPINV	common  inverse reduced gravity for bottom layer
!
!	HVM	selected horizontal velocity mode values
!	VVM	selected vertical velocity mode values
!	HCM	selected heat content mode values
!
!	AHMIX0	meqsf	depth of T=20 degrees
!	AHMIXP	meqsf	vertical velocity mode value at z=ahmix0
!

      CHARACTER ZIS*(*)
      REAL CPH(*),LIGHT(*),SLVL(*), VGRID(NZ)
      REAL HVM(NZ,*),VVM(NZ,*),HCM(NZ,*)
      REAL ahmix0(1),ahmixp(*),hvmsf(*)
      
!	EMHOR	horizontal velocity mode
!	VR	vertical velocity mode (not in common)
!	HC	heat content mode	
      REAL VR(MAXZVAL),EMHOR(MAXZVAL),HC(MAXZVAL)
      real THTEMP(1)
      
      REAL*8 VRDP(MAXZVAL) , WGHT , D1SQ,D2SQ,DZ2,SUMA
!	variable to pass grids to LININT
      REAL VIN(MAXZVAL)
      
!	decides on grid for interpolated variables
      IF(ZIS(1:1).EQ.'T')THEN
!	temperature grid
         DO 100 IZ = 1 , NPTS
            VIN(IZ) = TEMP1(IZ)
  100    CONTINUE
         
      ELSE
!	depth grid
         DO 105 IZ = 1 , NPTS
            VIN(IZ) = Z(IZ)
  105    CONTINUE
      ENDIF
!	output grid is based on ZTOP DZR NZ
      ZVAL = ZTOP
      
      
!     ******************************************************************
!     *  Calculates eigenvectors and prints out                        *
!     *   loops over modes                                             *
!     ******************************************************************
      
! CALCULATION OF HSTARP IN THE LOOP - 9.1.86 CHRIS
      HSTARP=0.
      
      DO 1400 JMODE = 2, NMODES
         IMODE = JMODE - 1
!	JMODE points to all the modes,
!		i.e., it keeps track of the data in the arrays.
!	IMODE points to the baraclinic modes, only.
         
         PHASEV = SQRT( ABS( 1.0 / EV(JMODE) ) )
         CPH(IMODE) = PHASEV
!
!  NEW CODE, 04 Feb 1986
         
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
!	WRITE(0,*)I,VRDP(I-1),NSQ(I),DZ(I),DZ(I+1)
  700    CONTINUE
         
!     Normalizes
!
!	Calculate horizontal mode
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
         
!	Normalization
!	If we write the vertical velocity eigenfunction as G(z),
!	then the horizontal velocity eigenfunction is Gz(z).
!	The normalization is such that
!	integral from -D to 0 of [Gz(z)] squared = D
!	In a constant buoyancy frequency ocean that choice means 
!	the maximum amplitude is the square root of 2.
!	The normalization also means that G(z) has units of meters.
!

         DO 995 J=1, NPTS
            EMHOR(J)=FAC*EMHOR(J)
  995    CONTINUE
         XDISP = PHASEV/NBAR
!	normalizes vertical velocity/displacement mode as well
         DO 996 J = 1 , NPTS
            VR(J) = FAC * VRDP(J)
  996    CONTINUE
         
!  CALCULATE LIGHTHILL COEFFICIENTS USING TRAPEZOIDAL RULE
!  SUM=INTEGRAL OF EMHOR FROM Z=0 TO Z=-DML
!      SUM = 0.0
!
!      DO 1310 J = 1, NPTS
!         IF ( Z(J+1) .GT. DML ) GO TO 1320
!         SUM = SUM + DZ(J+1) * (EMHOR(J+1)+EMHOR(J)) / 2.0
C1310  CONTINUE
!      GO TO 1330
!
C1320  SUM = SUM + EMHOR(J) * ( DML-Z(J) )
!      SUM = SUM + ( (EMHOR(J+1)-EMHOR(J)) * (DML-Z(J))**2 ) /
!     * ( 2. * DZ(J+1) )
!
C1330  SUMCL = SUM / DML
!	CL is the coupling coef. (M):
!	DU/DT + ... = TAUX/CL
!      CL = Z(NPTS) / ( SUMCL * EMHOR(1) )
!         
!	computes CL using just the surface value
         CL = Z(NPTS) / EMHOR(1)**2
         LIGHT(IMODE)=CL
!         print*,'lighthill coef (m)',imode,light(imode)
         HSTARP=HSTARP+1./CL
!	calculates SEALVL
         SLVL(IMODE) = EMHOR(1)/(EV(JMODE)*GRAV)
         
         HVMSF(IMODE) = EMHOR(1)
         
!	saves value at z0 = 50m (basic state profile depth)
! 	matches Steve's and Richard's TDEEP H input
         THTEMP(1) = 20.
!  really only need next line for one mode, oh well...
         CALL LININT (NPTS, TEMP1, Z, 1, THTEMP, AHMIX0, NLOW, NHIGH) 
!         print*,'calling LININT',AHMIX0
         CALL LININT (NPTS, TEMP1, VR, 1, THTEMP, AHMIXP(IMODE)
     *        ,NLOW, NHIGH)
!         print*,'calling LININT',VRDP(1),VRDP(2),VRDP(3),VRDP(4)
!         stop
         
         CALL LININT (NPTS , VIN , EMHOR ,NZ,VGRID, HVM(1,IMODE)
     *        ,NLOW , NHIGH )
         CALL LININT (NPTS , VIN , VR ,NZ,VGRID, VVM(1,IMODE)
     *        ,NLOW , NHIGH )

         
!	computes heat content
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
!	WRITE(0,*) DML,TLIGHT
!      print*,'TLIGHT = 1/HSTARP',TLIGHT

      return
      END

!      
!	linear interpolation to new grid
!	inputs
!	NT	number of points in the input profile
!	X	array of x
!	Y	array of y
!	NOUT
!	XOUT
!	outputs
!	NLOW	index of lowest XOUT within X
!	NHIGH	index of highest XOUT within X
!	YOUT	Y regridded to XOUT
!	
	subroutine LININT(NT,X,Y,NOUT,XOUT,YOUT,NLOW,NHIGH)
	REAL X(NT) , Y(NT) , XOUT(NOUT) , YOUT(NOUT)
!	initializes 
	XBINL = X(1)
	XBINH = X(NT)
	DELTAX = XBINH-XBINL
	DXOUT = XOUT(NOUT)-XOUT(1)
		IF((DXOUT*DELTAX).GT.0.0)THEN
!	both ascending or descending in X and XOUT
	IMIN = 1
	IMAX = NOUT
	ISTP = 1
				ELSE
!	opposite
	IMIN = NOUT
	IMAX = 1
	ISTP = -1
					ENDIF
	IFLIP = - ISTP

!	two versions -- one ascending in x, one descending
					IF(DELTAX.GT.0.0)THEN

	DO 20 I = IMIN , IMAX , ISTP
	IF(XBINL.LE.XOUT(I))GOTO 25
	YOUT(I) = Y(1)
20	CONTINUE
	I = 0
25	NLOW = I
!	finds NHIGH
	DO 30 I = IMAX , IMIN , IFLIP
	IF(XBINH.GE.XOUT(I))GOTO 35
	YOUT(I) = Y(NT)
30	CONTINUE
	I = 0
35	NHIGH = I
!	interpolates
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
!	descending
	DO 120 I = IMIN , IMAX , ISTP
	IF(XBINL.GE.XOUT(I))GOTO 125
	YOUT(I) = Y(1)
120	CONTINUE
	I = 0
125	NLOW = I
!	finds NHIGH
	DO 130 I = IMAX , IMIN , IFLIP
	IF(XBINH.LE.XOUT(I))GOTO 135
	YOUT(I) = Y(NT)
130	CONTINUE
	I = 0
135	NHIGH = I
!	interpolates
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
!	patch so that NLOW is always the lower index if the 
!		range is not empty
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
!
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,M,ISW,IER
      REAL*8               D(1001),E2(1001)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            II,I,JJ,J,K1,K
      REAL*8               DELTA,DLAM,EP,ERR,F,P,QP,Q,RLDLP,R,S,TOT
!                                  RLDLP = MACHINE PRECISION
      DATA               RLDLP/1.192093E-07/
!                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      DLAM = 0.0E0
      ERR = 0.0E0
      S = 0.0E0
!                                  LOOK FOR SMALL SUB-DIAGONAL ENTRIES
!                                  DEFINE INITIAL SHIFT FROM LOWER
!                                  GERSCHGORIN BOUND.
      TOT = D(1)
      Q = 0.0E0
      J = 0
      DO 15 I=1,N
         P = Q
         IF (I.EQ.1) GO TO 5
         IF (P.GT.RLDLP*(ABS(D(I))+ABS(D(I-1)))) GO TO 10
    5    E2(I) = 0.0E0
!                                  COUNT IF E2(I) HAS UNDERFLOWED
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
!                                  NEXT QR TRANSFORMATION
   35    TOT = TOT+S
         DELTA = D(N)-S
         I = N
         F = ABS(RLDLP*TOT)
         IF (DLAM.LT.F) DLAM = F
         IF (DELTA.GT.DLAM) GO TO 40
         IF (DELTA.GE.(-DLAM)) GO TO 75
         IER = 130
         GO TO 9000
!                                  REPLACE SMALL SUB-DIAGONAL SQUARES
!                                  BY ZERO TO REDUCE THE INCIDENCE OF
!                                  UNDERFLOWS
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
!                                  SET ERROR -- IRREGULAR END
!                                  DEFLATE MINIMUM DIAGONAL ELEMENT
         S = 0.0E0
         DELTA = QP
         DO 70 J=K,N
            IF (D(J).GT.DELTA) GO TO 70
            I = J
            DELTA = D(J)
   70    CONTINUE
!                                  CONVERGENCE
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
!   IMSL ROUTINE NAME   - UERTST
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - ELXSI/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
!
!   USAGE               - CALL UERTST (IER,NAME)
!
!   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
!                           IER = I+J WHERE
!                             I = 128 IMPLIES TERMINAL ERROR MESSAGE,
!                             I =  64 IMPLIES WARNING WITH FIX MESSAGE,
!                             I =  32 IMPLIES WARNING MESSAGE.
!                             J = ERROR CODE RELEVANT TO CALLING
!                                 ROUTINE.
!                NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING
!                           THE NAME OF THE CALLING ROUTINE. (INPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - UGETIO,USPKD
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
!                TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
!                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
!                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
!                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
!                UGETIO AS FOLLOWS..
!                                NIN = 0
!                                NOUT = NEW OUTPUT UNIT NUMBER
!                                CALL UGETIO(3,NIN,NOUT)
!                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      subroutine UERTST (IER,NAME)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      CHARACTER          NAME*(*)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEQDF,IOUNIT,LEVEL,LEVOLD,NIN,NMTB
      CHARACTER          IEQ,NAMEQ(6),NAMSET(6),NAMUPK(6)
      DATA               NAMSET/'U','E','R','S','E','T'/
      DATA               NAMEQ/6*' '/
      DATA               LEVEL/4/,IEQDF/0/,IEQ/'='/
!                                  UNPACK NAME INTO NAMUPK
!                                  FIRST EXECUTABLE STATEMENT
      CALL USPKD (NAME,6,NAMUPK,NMTB)
!                                  GET OUTPUT UNIT NUMBER
      CALL UGETIO(1,NIN,IOUNIT)
!                                  CHECK IER
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
!                                  PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
!                                  PRINT WARNING WITH FIX MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
      GO TO 30
   10 IF (IER.LE.32) GO TO 15
!                                  PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
      GO TO 30
   15 CONTINUE
!                                  CHECK FOR UERSET CALL
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
!                                  PRINT NON-DEFINED MESSAGE
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
!
!                                  SAVE P FOR P = R CASE
!                                    P IS THE PAGE NAMUPK
!                                    R IS THE ROUTINE NAMUPK
   55 IEQDF = 1
      DO 60 I=1,6
   60 NAMEQ(I) = NAMUPK(I)
   65 return
      END
!   IMSL ROUTINE NAME   - UGETIO
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - ELXSI/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1981
!
!   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
!                           VALUES FOR INPUT AND OUTPUT UNIT
!                           IDENTIFIERS.
!
!   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
!
!   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
!                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
!                           UNIT IDENTIFIER VALUES ARE returnED IN NIN
!                           AND NOUT, RESPECTIVELY.
!                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
!                           RESET FOR SUBSEQUENT USE.
!                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
!                           RESET FOR SUBSEQUENT USE.
!                NIN    - INPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
!                NOUT   - OUTPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
!                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
!                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
!                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
!                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      subroutine UGETIO(IOPT,NIN,NOUT)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/5/,NOUTD/6/
!                                  FIRST EXECUTABLE STATEMENT
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
!   IMSL ROUTINE NAME   - USPKD
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - ELXSI/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE
!                           CHARACTER STRING ARGUMENTS
!
!   USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
!
!   ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)
!                NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.
!                UNPAKD - CHARACTER ARRAY TO RECEIVE THE UNPACKED
!                         REPRESENTATION OF THE STRING. (OUTPUT)
!                NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - NONE
!
!   REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO A CHARACTER ARRAY
!                IN (A1) FORMAT.
!            2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF
!                THAT ARE IGNORED.
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
      subroutine USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,NCHARS,NCHMTB
!
      CHARACTER          UNPAKD(1),IBLANK
      CHARACTER*(1)      PACKED(1)
      DATA               IBLANK /' '/
!                                  INITIALIZE NCHMTB
      NCHMTB = 0
!                                  return IF NCHARS IS LE ZERO
      IF(NCHARS.LE.0) return
!                                  SET NC=NUMBER OF CHARS TO BE DECODED
      NC = MIN0 (129,NCHARS)
      DO 5 I=1,NC
         UNPAKD(I) = PACKED(I)
    5 CONTINUE
  150 FORMAT (129A1)
!                                  CHECK UNPAKD ARRAY AND SET NCHMTB
!                                  BASED ON TRAILING BLANKS FOUND
      DO 200 N = 1,NC
         NN = NC - N + 1
         IF(UNPAKD(NN) .NE. IBLANK) GO TO 210
  200 CONTINUE
  210 NCHMTB = NN
      return
      END
