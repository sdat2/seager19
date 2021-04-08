!-------------------------------------------------------------------------- 
      subroutine sst_update(nx,ny,delt,dxm,dym,hsfc,gamma,u1,v1,w1,t1,sstt,
     *            qnet,qflxt,dQdTt,dQdft,tdeep,udtdx,vdtdy,yv,idiag,qdiag,
     *            terms,ifper,mask)
      use wsst
      implicit none 
      include 'om_sst.h'
!-------------------------------------------------------------------------- 
!     SST advection equation on a A-GRID
!
!     delt     time step (s)
!     HSF!     depth of surface layer (m)
!     u1,v1    surface layer horizontal velocity components 
!     w1       surface layer vertical velocity component 
!     t1       model sst value
!     sstt     data sst value
!     qnet     net surface heat flux in watts/m**2
!     qflxt    data flux correction term
!     qdQdT    data dQ/dTs term
!     qdQdf    data dQ/df term
!     qdiag     diagnosed flux correction term
!     tdeep    temperature of deep water  
!     gamma    mixing fraction of sst and tdeep water
!     dxm      zonal grid spacing (m)
!     dym      merid grid spacing (m)
!
!    1. compute gradient of t1, with no flux boundary conditions
!
!    2. calculate w_up
!
!           w_up = max(0,w1)
!    where 
!             w1 = HSFC * div(u1,v1)
!
!    3. compute upwelling temperature
!
!        t_up = (1.-gamma)*t1 + gamma*tdeep
!
!    4. update sst
!
!  t1 = t1 + delt*(- u grad(t1) - w_up*(t1-t_up)/HSFC - qprime/HSFC 
!
!    where 
!         u grad(t1) = u1*d(t1)/dx + v1*d(t1)/dy 
!         qprime       = qnet/(Rho*Cp)
!         Cp  = 4184 watt-sec/(kg-degreeC)
!         rho = 1023 kg/m**3  (P=0,T=25,S=35)
!-------------------------------------------------------------------------- 
      integer nx, ny
      real  delt, hsfc, gamma, u1(nx,ny), v1(nx,ny), w1(nx,ny),
     *      tdeep(nx,ny), udtdx(nx,ny), vdtdy(nx,ny), dxm, dym,
     *      t1(nx,ny), qnet(nx,ny), sstt(nx,ny)
     *      ,yv(ny),qflxt(nx,ny),qdiag(nx,ny)
     *      ,terms(nx,ny,7)
     *      ,dQdTt(nx,ny),dQdft(nx,ny)

      integer i, j, k, ip, jp, im, jm, mask(nx,ny)
      logical ifper, idiag
      real RhoCp, uw, ue, ua, vs, vn, va, dudx, dvdy, dtdx, dtdy, td, sst,
     *     qprime, w_up, ugradt, qcon_inv, restore, flux_up

      RhoCp = 4184.*1023.   ! watt-sec/(m**3-degreeC)
      qcon_inv = 1./RhoCp

      do k = 1, npt_u
         i = ii(k) 
         j = jj(k) 
         ip = i + 1
         if (i .eq. nx) then
            if (ifper) then
               ip = 1
            else
               ip = i
            endif
         endif
         im = i - 1
         if (im .eq. 0) then
            if (ifper) then
               im = nx
            else
               im = i
            endif
         endif

         jp = min(ny,j + 1)
         jm = max(1 ,j - 1)

!    NOT STABLE, must do some upwinding
!         dtdx = (t1(ip,j) - t1(im,j))/2./dxm
!         dtdy = (t1(i,jp) - t1(i,jm))/2./dym

         if (u1(i,j) .gt. 0) then
            dtdx = (t1(i,j) - t1(im,j))/dxm
         else
            dtdx = (t1(ip,j) - t1(i,j))/dxm
         endif
         if (v1(i,j) .gt. 0) then
            dtdy = (t1(i,j) - t1(i,jm))/dym
         else
            dtdy = (t1(i,jp) - t1(i,j))/dym
         endif

         udtdx(i,j) = u1(i,j) * dtdx
         vdtdy(i,j) = v1(i,j) * dtdy

      enddo

      do i = 1, nx
      do j = 1, ny
         if (mask(i,j) .eq. 0) then
            t1(i,j) = sstt(i,j)
            qnet(i,j) = 0.
         else
            sst = t1(i,j)

            td = MIN(tdeep(i,j),sst-.5)
            w_up  = max(0.,w1(i,j))
            
            ugradt = udtdx(i,j) + vdtdy(i,j)
            flux_up = w_up* gamma* (sst-td)/ HSFC
            qprime = qcon_inv* qnet(i,j)

            terms(i,j,1) = -udtdx(i,j)
            terms(i,j,2) = -vdtdy(i,j)
            terms(i,j,3) = -flux_up
            terms(i,j,4) = -qprime
            terms(i,j,5) = qflxt(i,j)/HSFC
            terms(i,j,6) = td
            terms(i,j,7) = w_up
            t1(i,j) = sst - delt*(ugradt + flux_up + qprime) 
     *                    + delt*qflxt(i,j)/HSFC 

!  next two lines only for diagnosing the qflx
            if (idiag) then
               qdiag(i,j) = HSFC*(sstt(i,j) - t1(i,j))/delt
               t1(i,j) = sstt(i,j)
            endif

         endif
         
      enddo 
      enddo 

      return
      end

!-------------------------------------------------------------------------- 
      subroutine tdeep_calc(nx,ny,nmodes,ahmix0,ahmixp,hm,tdeep,hmod,htherm)
      use wsst
      implicit none 
      include 'om_sst.h'
!-------------------------------------------------------------------------- 
!      Subroutine for calculating tdeep (temperature of upwelled
!      water) given thermocline depths
!-------------------------------------------------------------------------- 
      integer nx,ny,nmodes
      real ahmixp(nmodes,nx,ny), hm(ny,nx,nmodes), tdeep(nx,ny)
      real hmod(nx,ny), htherm(nx,ny), ahmix0(nx,ny)

      integer NKNOTSMAX
      PARAMETER(NKNOTSMAX=20)
      REAL*8 spknots(NKNOTSMAX), spcoeffs(NKNOTSMAX)
      REAL B0,B1,B2,HCUT,TCUT
      INTEGER nknots,nspk
      COMMON /SSTPARAM/ spknots,spcoeffs,B0,B1,B2,HCUT,TCUT,nspk,nknots

      real hmodel, tdeeph, h20
      integer i,j,k,im
      real tdeep0, tdeep1, h20_0, h20_1

      do k = 1, npt_u
         i = ii(k)
         j = jj(k)
         hmodel = ahmix0(i,j)
         do im= 1, nmodes
            hmodel= hmodel + ahmixp(im,i,j)*hm(j,i,im)
         enddo
      
! looks like Benno did fit on (hmodel-ahmix0) - I have changed this
         h20=b0 + b1*hmodel + b2*hmodel**2
      
!  shouldn't be needed, but for global domain the peruvian coast
!     upwelling is causing trouble:
         h20 = max(h20,20.0)

         hmod(i,j) = hmodel
         htherm(i,j) = h20
         tdeep(i,j) = tdeeph(hmodel)

!         TDEEP = MIN(TDEEP,SST-.5)
! old way was to fit hmodel to h20 and then h20 to tdeep
! orig    tdeep(i,j) = tdeeph_ERA(h20)  
!
!  could use a linear fit instead:
!         h20_0  = 60
!         tdeep0 = 24
!         h20_1  = 180
!         tdeep1 = 28.4
!         tdeep(i,j) = tdeep0 + hmodel* (tdeep1- tdeep0)/ (h20_1- h20_0)
!
! we're not using this for now:
!
      enddo

      return
      end


!-------------------------------------------------------------------------- 
      real function tdeeph(H20)
!-------------------------------------------------------------------------- 

      integer NKNOTSMAX
      PARAMETER(NKNOTSMAX=20)
      REAL*8 spknots(NKNOTSMAX), spcoeffs(NKNOTSMAX)
      REAL B0,B1,B2,HCUT,TCUT
      INTEGER nknots,nspk
      COMMON /SSTPARAM/ spknots,spcoeffs,B0,B1,B2,HCUT,TCUT,nspk,nknots

      REAL*4 H20
      REAL*8 x(1), y(1)
      integer e

      e = 0
      x(1) = H20

      IF(H20.LT.HCUT)THEN
         tdeeph=TCUT
      ELSE
         call splev(spknots(1:nknots),nknots,spcoeffs(1:nknots),nspk,x,y,1,e,ier)
         tdeeph = y(1)
      ENDIF

      return
      end

!-------------------------------------------------------------------------- 
      subroutine sst_mem (nxy) 
      use wsst
      implicit none  
      include 'om_sst.h' 
!-------------------------------------------------------------------------- 
      integer nxy

      allocate (ii(nxy))
      allocate (jj(nxy))

      return
      end

!-------------------------------------------------------------------------- 
      subroutine sst_init(nx,ny,sst,sstt,xv,yv,tp,irest)
      use wsst
      implicit none 
      include 'om_sst.h'
!-------------------------------------------------------------------------- 
      integer nx, ny
      real sst(nx,ny), sstt(nx,ny),tp(nx,ny)
      real xv(nx),yv(ny),cy,cx

      integer i, j, k, irest

      if (irest.le.2) then
         do k = 1, npt_u
            i = ii(k)
            j = jj(k)
            sst(i,j) = sstt(i,j)
         enddo
      endif

      return
      end
      
!-------------------------------------------------------------------------- 
      subroutine sst_index(nx,ny)
      use wsst
      implicit none 
      include 'om_core.h'  ! for npt,ival,jval
      include 'om_sst.h'
!-------------------------------------------------------------------------- 
      integer nx, ny

      integer  k, nxy

      nxy = nx*ny
      npt_u = npt

      call sst_mem(nxy)

      do k = 1, npt_u
         ii(k) = ival(k)
         jj(k) = jval(k)
      enddo

      return
      end
      
!-------------------------------------------------------------------------- 
      subroutine upwelling(nx,ny,dxm,dym,hsfc,u1,v1,w1,ifper)
      use wsst
      implicit none 
      include 'om_sst.h'
!-------------------------------------------------------------------------- 
!     HSF!     depth of surface layer (m)
!     u1,v1    surface layer horizontal velocity components 
!     w1       surface layer vertical velocity component (on scalar points)
!     dxm      zonal grid spacing (m)
!     dym      merid grid spacing (m)
!
!    compute upwelling velocity
!
!             w1 = HSFC * div(u1,v1)
!
!-------------------------------------------------------------------------- 
      integer nx, ny
      real  hsfc, u1(nx,ny), v1(nx,ny), w1(nx,ny),dxm, dym
      logical ifper

      integer i, j, k, ip, jp, im, jm
      real uw, ue, vs, vn, dudx, dvdy

      do k = 1, npt_u
         i = ii(k) 
         j = jj(k) 
         ip = i+1
         if (i .eq. nx) then
            if (ifper) then
               ip = 1
            else
               ip = i
            endif
         endif
         im = i - 1
         if (im .eq. 0) then
            if (ifper) then
               im = nx
            else
               im = i
            endif
         endif

         jp = min(ny,j + 1)
         jm = max(1, j - 1)

         uw = (u1(im,j) + u1(i,j))/2.
         ue = (u1(ip,j) + u1(i,j))/2.
         dudx = (ue - uw)/dxm

         vs = (v1(i,jm) + v1(i,j))/2.
         vn = (v1(i,jp) + v1(i,j))/2.
         dvdy = (vn - vs)/dym
    
         w1(i,j) = HSFC*(dudx + dvdy)
      enddo
    
      return
      end


      subroutine splev(t,n,c,k,x,y,m,e,ier)
!  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
!  a spline s(x) of degree k, given in its b-spline representation.
!
!  calling sequence:
!     call splev(t,n,c,k,x,y,m,e,ier)
!
!  input parameters:
!    t    : array,length n, which contains the position of the knots.
!    n    : integer, giving the total number of knots of s(x).
!    !    : array,length n, which contains the b-spline coefficients.
!    k    : integer, giving the degree of s(x).
!    x    : array,length m, which contains the points where s(x) must
!           be evaluated.
!    m    : integer, giving the number of points where s(x) must be
!           evaluated.
!    e    : integer, if 0 the spline is extrapolated from the end
!           spans for points not in the support, if 1 the spline
!           evaluates to zero for those points, if 2 ier is set to
!           1 and the subroutine returns, and if 3 the spline evaluates
!           to the value of the nearest boundary point.
!
!  output parameter:
!    y    : array,length m, giving the value of s(x) at the different
!           points.
!    ier  : error flag
!      ier = 0 : normal return
!      ier = 1 : argument out of bounds and e == 2
!      ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!    m >= 1
!--    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl.
!
!  references :
!    de boor c  : on calculating with b-splines, j. approximation theory
!                 6 (1972) 50-62.
!    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
!                 applics 10 (1972) 134-149.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!++ pearu: 11 aug 2003
!++   - disabled cliping x values to interval [min(t),max(t)]
!++   - removed the restriction of the orderness of x values
!++   - fixed initialization of sp to double precision value
!
!  ..scalar arguments..
      integer n, k, m, e, ier
!  ..array arguments..
      real*8 t(n), c(n), x(m), y(m)
!  ..local scalars..
      integer i, j, k1, l, ll, l1, nk1
!++..
      integer k2
!..++
      real*8 arg, sp, tb, te
!  ..local array..
      real*8 h(20)
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
!--      if(m-1) 100,30,10
!++..
      if (m .lt. 1) go to 100
!..++
!--  10  do 20 i=2,m
!--        if(x(i).lt.x(i-1)) go to 100
!--  20  continue
      ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
      k1 = k + 1
!++..
      k2 = k1 + 1
!..++
      nk1 = n - k1
      tb = t(k1)
      te = t(nk1 + 1)
      l = k1
      l1 = l + 1
!  main loop for the different points.
      do 80 i = 1, m
!  fetch a new x-value arg.
        arg = x(i)
!  check if arg is in the support
        if (arg .lt. tb .or. arg .gt. te) then
            if (e .eq. 0) then
                goto 35
            else if (e .eq. 1) then
                y(i) = 0
                goto 80
            else if (e .eq. 2) then
                ier = 1
                goto 100
            else if (e .eq. 3) then
                if (arg .lt. tb) then
                    arg = tb
                else
                    arg = te
                endif
            endif
        endif
!  search for knot interval t(l) <= arg < t(l+1)
!++..
 35     if (arg .ge. t(l) .or. l1 .eq. k2) go to 40
        l1 = l
        l = l - 1
        go to 35
!..++
  40    if(arg .lt. t(l1) .or. l .eq. nk1) go to 50
        l = l1
        l1 = l + 1
        go to 40
!  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t, n, k, arg, l, h)
!  find the value of s(x) at x=arg.
        sp = 0.0d0
        ll = l - k1
        do 60 j = 1, k1
          ll = ll + 1
          sp = sp + c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end

      subroutine fpbspl(t,n,k,x,l,h)
!  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
!  degree k at t(l) <= x < t(l+1) using the stable recurrence
!  relation of de boor and cox.
!  Travis Oliphant  2007
!    changed so that weighting of 0 is used when knots with
!      multiplicity are present.
!    Also, notice that l+k <= n and 1 <= l+1-k
!      or else the routine will be accessing memory outside t
!      Thus it is imperative that that k <= l <= n-k but this
!      is not checked.
!  ..
!  ..scalar arguments..
      real*8 x
      integer n,k,l
!  ..array arguments..
      real*8 t(n),h(20)
!  ..local scalars..
      real*8 f,one
      integer i,j,li,lj
!  ..local arrays..
      real*8 hh(19)
!  ..
      one = 0.1d+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.0d0
        do 20 i=1,j
          li = l+i
          lj = li-j
          if (t(li).ne.t(lj)) goto 15
          h(i+1) = 0.0d0
          goto 20
  15      f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end

