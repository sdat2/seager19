c-------------------------------------------------------------------------- 
      subroutine sst_update(nx,ny,delt,dxm,dym,hsfc,gamma,u1,v1,w1,t1,sstt,
     *            qnet,qflxt,dQdTt,dQdft,tdeep,udtdx,vdtdy,yv,idiag,qdiag,
     *            terms,ifper,mask)
      use wsst
      implicit none 
      include 'om_sst.h'
c-------------------------------------------------------------------------- 
c     SST advection equation on a A-GRID
c
c     delt     time step (s)
c     HSFC     depth of surface layer (m)
c     u1,v1    surface layer horizontal velocity components 
c     w1       surface layer vertical velocity component 
c     t1       model sst value
c     sstt     data sst value
c     qnet     net surface heat flux in watts/m**2
c     qflxt    data flux correction term
c     qdQdT    data dQ/dTs term
c     qdQdf    data dQ/df term
c     qdiag     diagnosed flux correction term
c     tdeep    temperature of deep water  
c     gamma    mixing fraction of sst and tdeep water
c     dxm      zonal grid spacing (m)
c     dym      merid grid spacing (m)
c
c    1. compute gradient of t1, with no flux boundary conditions
c
c    2. calculate w_up
c
c           w_up = max(0,w1)
c    where 
c             w1 = HSFC * div(u1,v1)
c
c    3. compute upwelling temperature
c
c        t_up = (1.-gamma)*t1 + gamma*tdeep
c
c    4. update sst
c
c  t1 = t1 + delt*(- u grad(t1) - w_up*(t1-t_up)/HSFC - qprime/HSFC 
c
c    where 
c         u grad(t1) = u1*d(t1)/dx + v1*d(t1)/dy 
c         qprime       = qnet/(Rho*Cp)
c         Cp  = 4184 watt-sec/(kg-degreeC)
c         rho = 1023 kg/m**3  (P=0,T=25,S=35)
c-------------------------------------------------------------------------- 
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

c    NOT STABLE, must do some upwinding
c         dtdx = (t1(ip,j) - t1(im,j))/2./dxm
c         dtdy = (t1(i,jp) - t1(i,jm))/2./dym

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

c  next two lines only for diagnosing the qflx
            if (idiag) then
               qdiag(i,j) = HSFC*(sstt(i,j) - t1(i,j))/delt
               t1(i,j) = sstt(i,j)
            endif

         endif
         
      enddo 
      enddo 

      return
      end

c-------------------------------------------------------------------------- 
      subroutine tdeep_calc(nx,ny,nmodes,ahmix0,ahmixp,hm,tdeep,hmod,htherm)
      use wsst
      implicit none 
      include 'om_sst.h'
c-------------------------------------------------------------------------- 
c      Subroutine for calculating tdeep (temperature of upwelled
c      water) given thermocline depths
c-------------------------------------------------------------------------- 
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
      
c looks like Benno did fit on (hmodel-ahmix0) - I have changed this
         h20=b0 + b1*hmodel + b2*hmodel**2
      
c  shouldn't be needed, but for global domain the peruvian coast
c     upwelling is causing trouble:
         h20 = max(h20,20.0)

         hmod(i,j) = hmodel
         htherm(i,j) = h20
         tdeep(i,j) = tdeeph(hmodel)

c         TDEEP = MIN(TDEEP,SST-.5)
c old way was to fit hmodel to h20 and then h20 to tdeep
c orig    tdeep(i,j) = tdeeph_ERA(h20)  
c
c  could use a linear fit instead:
c         h20_0  = 60
c         tdeep0 = 24
c         h20_1  = 180
c         tdeep1 = 28.4
c         tdeep(i,j) = tdeep0 + hmodel* (tdeep1- tdeep0)/ (h20_1- h20_0)
c
c we're not using this for now:
c
      enddo

      return
      end


c-------------------------------------------------------------------------- 
      real function tdeeph(H20)
c-------------------------------------------------------------------------- 

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

c-------------------------------------------------------------------------- 
      subroutine sst_mem (nxy) 
      use wsst
      implicit none  
      include 'om_sst.h' 
c-------------------------------------------------------------------------- 
      integer nxy

      allocate (ii(nxy))
      allocate (jj(nxy))

      return
      end

c-------------------------------------------------------------------------- 
      subroutine sst_init(nx,ny,sst,sstt,xv,yv,tp,irest)
      use wsst
      implicit none 
      include 'om_sst.h'
c-------------------------------------------------------------------------- 
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
      
c-------------------------------------------------------------------------- 
      subroutine sst_index(nx,ny)
      use wsst
      implicit none 
      include 'om_core.h'  ! for npt,ival,jval
      include 'om_sst.h'
c-------------------------------------------------------------------------- 
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
      
c-------------------------------------------------------------------------- 
      subroutine upwelling(nx,ny,dxm,dym,hsfc,u1,v1,w1,ifper)
      use wsst
      implicit none 
      include 'om_sst.h'
c-------------------------------------------------------------------------- 
c     HSFC     depth of surface layer (m)
c     u1,v1    surface layer horizontal velocity components 
c     w1       surface layer vertical velocity component (on scalar points)
c     dxm      zonal grid spacing (m)
c     dym      merid grid spacing (m)
c
c    compute upwelling velocity
c
c             w1 = HSFC * div(u1,v1)
c
c-------------------------------------------------------------------------- 
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
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,e,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c    e    : integer, if 0 the spline is extrapolated from the end
c           spans for points not in the support, if 1 the spline
c           evaluates to zero for those points, if 2 ier is set to
c           1 and the subroutine returns, and if 3 the spline evaluates
c           to the value of the nearest boundary point.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier = 1 : argument out of bounds and e == 2
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c--    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c++ pearu: 11 aug 2003
c++   - disabled cliping x values to interval [min(t),max(t)]
c++   - removed the restriction of the orderness of x values
c++   - fixed initialization of sp to double precision value
c
c  ..scalar arguments..
      integer n, k, m, e, ier
c  ..array arguments..
      real*8 t(n), c(n), x(m), y(m)
c  ..local scalars..
      integer i, j, k1, l, ll, l1, nk1
c++..
      integer k2
c..++
      real*8 arg, sp, tb, te
c  ..local array..
      real*8 h(20)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
c--      if(m-1) 100,30,10
c++..
      if (m .lt. 1) go to 100
c..++
c--  10  do 20 i=2,m
c--        if(x(i).lt.x(i-1)) go to 100
c--  20  continue
      ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k + 1
c++..
      k2 = k1 + 1
c..++
      nk1 = n - k1
      tb = t(k1)
      te = t(nk1 + 1)
      l = k1
      l1 = l + 1
c  main loop for the different points.
      do 80 i = 1, m
c  fetch a new x-value arg.
        arg = x(i)
c  check if arg is in the support
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
c  search for knot interval t(l) <= arg < t(l+1)
c++..
 35     if (arg .ge. t(l) .or. l1 .eq. k2) go to 40
        l1 = l
        l = l - 1
        go to 35
c..++
  40    if(arg .lt. t(l1) .or. l .eq. nk1) go to 50
        l = l1
        l1 = l + 1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t, n, k, arg, l, h)
c  find the value of s(x) at x=arg.
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
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  Travis Oliphant  2007
c    changed so that weighting of 0 is used when knots with
c      multiplicity are present.
c    Also, notice that l+k <= n and 1 <= l+1-k
c      or else the routine will be accessing memory outside t
c      Thus it is imperative that that k <= l <= n-k but this
c      is not checked.
c  ..
c  ..scalar arguments..
      real*8 x
      integer n,k,l
c  ..array arguments..
      real*8 t(n),h(20)
c  ..local scalars..
      real*8 f,one
      integer i,j,li,lj
c  ..local arrays..
      real*8 hh(19)
c  ..
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

