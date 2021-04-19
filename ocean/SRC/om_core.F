c----------------------------------------------------------------------
      subroutine inc_init (nxr,nyr,nxyc,k,mask,dt,hx,hy,fc,rayl,rlap,
     *                     um,vm,hm,un,vn,hn,vmy,ischeme,ifper)

      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
c     this subroutine precomputes various quantities, indexing stuff
c          and coefficients for tridiagonal solves

c     input:
c        k       : vertical mode number
c        nxr     : number of grid points in zonal direction
c        nyr     : number of grid points in meridional direction
c        nxyc    : number of water points as calculated by make_iox
c        mask    : land mask (0=land, nonzero=ocean, or edge of ocean)
c        dt      : time step in scaled time units
c        hx      : zonal grid spacing (scaled)
c        hy      : meridional grid spacing (scaled)
c        fc      : scaled coriolis parameter
c        rayl    : scaled Rayleigh coefficients
c        rlap    : scaled Laplacian diffusion coefficients
c        hr      : initial pressure  (scaled)
c        ur      : initial zonal velocity (scaled)
c        vr      : initial meridional velocity (scaled)
c        txr     : initial zonal momentum forcing term (scaled)
c        tyr     : initial meridional momentum forcing term (scaled)
      common /constants/ xconv,yconv,rl

      real xconv,yconv,delx,dely,rl
      logical ifper
      integer nxr,nyr,k,nxyc,ischeme
      integer mask(nxr,nyr)
      real fc(nyr), dt, hx, hy, rayl, rlap
      real um(nyr,nxr),vm(nyr,nxr),vmy(nyr,nxr),hm(nyr,nxr)
      real un(nyr,nxr),vn(nyr,nxr),hn(nyr,nxr)

      integer i,j,ip,jp
      real fcc,fp,fm
      real d1,d2,db,dp
      integer in,ipt
      real tausqi,csp,bb,xc,yc,xl,yl

      tau(k) = dt/2.0
      dt2(k) = 2.0*dt
      taui(k) = 2.0/dt
      tausqi =  taui(k)**2
      thxi(k) = 1.0 / ( hx + hx)
      thyi(k) = 1.0 / ( hy + hy)
      hxsq(k) = hx*hx
      hysq(k) = hy*hy
      eps(k) = rayl*tau(k)
      vnu(k) = rlap

      call make_indis(nxr,nyr,nxyc,k,mask,dt,hx,ifper,ischeme)

      do j = 1, nyr
         fcor(j,k) = fc(j)
         do i = 1, nxr
            gamm(j,i,k)=  .0
            gamp(j,i,k)=  .0
            gamc(j,i,k)= 1.0
         enddo
      enddo
      fcor(nyr+1,k) = 2.*fcor(nyr,k) - fcor(nyr-1,k)
      
      do in = 1 , nptix
         ipt = lptix(in)
         i   = ival(ipt)
         j   = jval(ipt)
         fcc = fcor(j,k)
         d1 = dt_l(i,j,k)
         d2 = dt_r(i,j,k)
         db = tau(k)*(d1+d2)/2.0
         dp = tau(k)*(d1-d2)/2.0
         gamm(j,i,k)  =  1.0
         gamp(j,i,k)  =  1.0
         gamc(j,i,k)  = -2.0 
     *            - hysq(k)*(fcc**2 + taui(k)/db + dp/db)
      enddo
      
      if (ischeme.eq.1) then
         do in = 1 , npt_w
            ipt = lpt_w(in)
            i   = ival(ipt)
            j   = jval(ipt)
            fcc      = fcor(j,k)
            fp       = fcor(j+1,k)
            fm       = fcor(j-1,k)
            gamm(j,i,k)=  1.0 - hysq(k)*(fm - abs(fm))*thyi(k)
            gamp(j,i,k)=  1.0 + hysq(k)*(fp + abs(fp))*thyi(k)
            gamc(j,i,k)= -2.0 - hysq(k)*(tausqi + 2.*abs(fcc)*thyi(k))
         enddo
         
         do in = 1 , npt_e
            ipt = lpt_e(in)
            i   = ival(ipt)
            j   = jval(ipt)
            fcc      = fcor(j,k)
            fp       = fcor(j+1,k)
            fm       = fcor(j-1,k)
            gamm(j,i,k)=  1.0 + hysq(k)*(fm + abs(fm))*thyi(k)
            gamp(j,i,k)=  1.0 - hysq(k)*(fp - abs(fp))*thyi(k)
            gamc(j,i,k)= -2.0 - hysq(k)*(tausqi + 2.*abs(fcc)*thyi(k))
         enddo

      else
         do in = 1 , npt_w
            ipt = lpt_w(in)
            i   = ival(ipt)
            j   = jval(ipt)
            gamm(j,i,k)  =  1.0 
            gamp(j,i,k)  =  1.0 
            gamc(j,i,k)  = -2.0 - hysq(k)*tausqi
         enddo
         
         do in = 1 , npt_e
            ipt = lpt_e(in)
            i   = ival(ipt)
            j   = jval(ipt)
            gamm(j,i,k)  =  1.0
            gamp(j,i,k)  =  1.0
            gamc(j,i,k)  = -2.0 - hysq(k)*tausqi
         enddo
      endif
      
      do in = 1 , npt_n
         ipt = lpt_n(in)
         i   = ival(ipt)
         j   = jval(ipt)
         gamm(j,i,k)  = 0.0
         gamp(j,i,k)  = 0.0
         gamc(j,i,k)  = 1.0
      enddo
      
      do in = 1 , npt_s
         ipt = lpt_s(in)
         i   = ival(ipt)
         j   = jval(ipt)
         gamm(j,i,k)  = 0.0
         gamp(j,i,k)  = 0.0
         gamc(j,i,k)  = 1.0
      enddo

      call bcset(nxr,nyr,um,vm)
      call make_vy(1,thyi(k),nxr,nyr,vm,vmy)

c  set up for leapfrog

      if (ischeme.eq.2) then
         do ipt = 1 , npt
            i   = ival(ipt)
            j   = jval(ipt)
            hn(j,i) = hm(j,i)
            un(j,i) = um(j,i)
            vn(j,i) = vm(j,i)
         enddo
      endif

      return
      end
      
c----------------------------------------------------------------------
      subroutine inc (nstep,k,nxr,nyr,hm,um,vm,vmy,F,G,Q,ifper)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
c
c     solves:
c
c        d(um)/dt - fcor*vm  + d(hm)/dx = F - r*um + nu*(u_xx + u_yy)
c
c        d(vm)/dt + fcor*um  + d(hm)/dy = G - r*vm + nu*(v_xx + v_yy)
c
c        d(hm)/dt + d(um)/dx + d(vm)/dy = Q - r*hm + nu*(h_xx + h_yy)
c
c      on an A-grid


c     input:
c        nxr      : number of grid points in zonal direction
c        nyr      : number of grid points in meridional direction
c        hm,um,vm : k-th scaled mode
c        F,G      : momentum forcing
c        Q        : mass forcing

      logical ifper
      integer nstep, k, nxr, nyr
      real F(nyr,nxr),G(nyr,nxr),Q(nyr,nxr)
      real hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr),vmy(nyr,nxr)

      integer i,j,ipt

c  add terms at the old time
      call add_diff(k,nxr,nyr,hm,um,vm,F,G,Q)

      call make_rhs4uh(k,nxr,nyr,hm,um,vm,vmy,F,Q)
      call add2rhs4uh(k,nxr,nyr,F,Q)

c  add terms at the new time
      call make_rhs4v(k,nxr,nyr,hm,um,vm,G)
      call add2rhs4v(k,nxr,nyr,G)
      
      call solve_v(nxr,nyr,NXMAX,NYMAX, rhs, vm
     *            , gamm(1,1,k), gamc(1,1,k), gamp(1,1,k))

      call make_vy(1,thyi(k),nxr,nyr,vm,vmy)
      call solve_uh(k,nxr,nyr,hm,um,vm,vmy)

      if (ifper) then
         do j = 1, nyr
            hm(j,1)   = hm(j,nxr-3)
            hm(j,2)   = hm(j,nxr-2)
            hm(j,nxr  ) = hm(j,4)
            hm(j,nxr-1) = hm(j,3)
            um(j,1)   = um(j,nxr-3)
            um(j,2)   = um(j,nxr-2)
            um(j,nxr  ) = um(j,4)
            um(j,nxr-1) = um(j,3)
            vm(j,1)   = vm(j,nxr-3)
            vm(j,2)   = vm(j,nxr-2)
            vm(j,nxr  ) = vm(j,4)
            vm(j,nxr-1) = vm(j,3)
         enddo
      endif

      call filt(nstep,nxr,nyr,hm,um,vm,vmy,2.*eps(k))

      if (ifper) then
         do j = 1, nyr
            hm(j,1)   = hm(j,nxr-3)
            hm(j,2)   = hm(j,nxr-2)
            hm(j,nxr  ) = hm(j,4)
            hm(j,nxr-1) = hm(j,3)
            um(j,1)   = um(j,nxr-3)
            um(j,2)   = um(j,nxr-2)
            um(j,nxr  ) = um(j,4)
            um(j,nxr-1) = um(j,3)
            vm(j,1)   = vm(j,nxr-3)
            vm(j,2)   = vm(j,nxr-2)
            vm(j,nxr  ) = vm(j,4)
            vm(j,nxr-1) = vm(j,3)
            vmy(j,1)   = vmy(j,nxr-3)
            vmy(j,2)   = vmy(j,nxr-2)
            vmy(j,nxr  ) = vmy(j,4)
            vmy(j,nxr-1) = vmy(j,3)
         enddo
      endif


      return
      end
      
c----------------------------------------------------------------------
      subroutine make_indis(nx,ny,nxyc,k,mask,dt,hx,ifper,ischeme)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nx,ny,mask(nx,ny),k,nxyc,ischeme
      real dt,hx

      real shapi_coef, shapb_coef, shaph_coef, shapi_coefx, shapi_coefy
      integer nshapxd, nshapyd, nshapb, nshap, nshapt
      common /filter/ nshap,nshapb,nshapxd,nshapyd,shapi_coef,shapb_coef
     *               , shaph_coef, nshapt, shapi_coefx, shapi_coefy

c    NB: - mask is assumed to be defined on an A-grid which is the one
c           in LOAM (in a cell-centered view, the land intersects
c                 the cells through the middle; creating half and 
c                 quarter sized cells around the boundary of the ocean)
c           therefore, each water segment must be at least 2 grid points wide

      logical ifper
      integer i,j,mw,me,ms,mn,idtx,ipt,ix,is
      integer mm,m0,i0,i1,i2,i3,m3,ind,ndx,ndy,in,idpt,id,jd
      real alf

      if (nx.gt.NXMAX.or.ny.gt.NYMAX) then
         print*,'nx or ny is too large, redefine NXMAX or NYMAX in meqsf.h'
         stop
      endif

      npt = 0
      do j = 1, ny
         do i = 1, nx
            dt_l(i,j,k) = 1.
            dt_r(i,j,k) = 1.
            it_l(i,j,k) = i
            it_r(i,j,k) = i
            at_l(i,j,k) = 0.
            at_r(i,j,k) = 0.
            if (mask(i,j).ne.0) then
               npt = npt + 1
               ival(npt) = i
               jval(npt) = j
            endif
         enddo
      enddo

      if (npt.ne.nxyc) then
         print*,'npt.ne.nxyc in make_indis'
         stop
      endif

      npt_w = 0
      npt_e = 0
      npt_ew = 0
      npt_n = 0
      npt_s = 0
      npt_ns = 0
      nptix = 0
      nptiy = 0
      nptixy = 0

      do ipt = 1, npt
         i = ival(ipt)
         j = jval(ipt)
         mw = 0
         me = 0
         ms = 0
         mn = 0
         if (i.gt.1)  mw = mask(i-1,j)
         if (i.lt.nx) me = mask(i+1,j)
         if (j.gt.1)  ms = mask(i,j-1)
         if (j.lt.ny) mn = mask(i,j+1)
            
         if (mw.eq.0) then
c                add grid point to list of western boundary points
            npt_w = npt_w + 1
            lpt_w(npt_w) = ipt
c               print*,'pt_w',npt_w,ipt
         endif
         if (me.eq.0) then
c                add grid point to list of eastern boundary points
            npt_e = npt_e + 1
            lpt_e(npt_e) = ipt
         endif
         if (mn.eq.0) then
c                add grid point to list of northern boundary points
            npt_n = npt_n + 1
            lpt_n(npt_n) = ipt
         endif
         if (ms.eq.0) then
c                add grid point to list of southern boundary points
            npt_s = npt_s + 1
            lpt_s(npt_s) = ipt
         endif
         if (mw*me.ne.0) then
c                add grid point to list of interior (x) points
            nptix = nptix + 1
            lptix(nptix) = ipt
         else
c                add grid point to list of boundary (x) points
            npt_ew = npt_ew + 1
            lpt_ew(npt_ew) = ipt
         endif
         if (ms*mn.ne.0) then
c                add grid point to list of interior (y) points
            nptiy = nptiy + 1
            lptiy(nptiy) = ipt
         else
c                add grid point to list of boundary (y) points
            npt_ns = npt_ns + 1
            lpt_ns(npt_ns) = ipt
         endif
         if (ms*mn*me*mw.ne.0) then
c                add grid point to list of interior points
            nptixy = nptixy + 1
            lptixy(nptixy) = ipt
         else
c                add grid point to list of boundary points
            npt_b = npt_b + 1
            lpt_b(npt_b) = ipt
            if (ms+me.eq.0.or.ms+mw.eq.0.or.mn+me.eq.0.or.mn+mw.eq.0) then
               npt_c = npt_c + 1
               lpt_c(npt_c) = ipt
            else
               npt_bc = npt_bc + 1
               lpt_bc(npt_bc) = ipt
            endif
         endif
         if (me*mw.ne.0) then
            if (mn.eq.0) then
               npt_no = npt_no + 1
               lpt_no(npt_no) = ipt
            elseif (ms.eq.0) then
               npt_so = npt_so + 1
               lpt_so(npt_so) = ipt
            endif
         endif
         if (ms*mn.ne.0) then
            if (me.eq.0) then
               npt_eo = npt_eo + 1
               lpt_eo(npt_eo) = ipt
            elseif (mw.eq.0) then
               npt_wo = npt_wo + 1
               lpt_wo(npt_wo) = ipt
            endif
         endif
         if (mn+mw.eq.0) then
            npt_nwo = npt_nwo + 1
            lpt_nwo(npt_nwo) = ipt
         endif
         if (mn+me.eq.0) then
            npt_neo = npt_neo + 1
            lpt_neo(npt_neo) = ipt
         endif
         if (ms+mw.eq.0) then
            npt_swo = npt_swo + 1
            lpt_swo(npt_swo) = ipt
         endif
         if (ms+me.eq.0) then
            npt_seo = npt_seo + 1
            lpt_seo(npt_seo) = ipt
         endif
      enddo

      idtx = int(dt/hx)
      alf = dt/hx - idtx

c     find params for left characteristics of interior points
      nrt_lc(k) = 0
      nbt_lc(k) = 0
      not_lc(k) = 0
      do ix = 1, nptix
         ipt = lptix(ix)
         i   = ival(ipt)
         j   = jval(ipt)
         do is = 1 , idtx
            mm = 0
            if (i-is-1.ge.1) mm = mask(i-is-1,j)
            if (mm.eq.0) then  !characteristic intersects boundary
               not_lc(k) = not_lc(k) + 1
               lot_lc(not_lc(k),k) = ipt
               it_l(i,j,k) = i - is  ! i value of boundary point
               dt_l(i,j,k) = is*hx/dt
               goto 10
            endif
         enddo
c            characteristic did not intersect boundary
c               takeoff point is in the interval:  (i1,i2]
         i2 = i - idtx
         i1 = i2 - 1
         i0 = i2 - 2
         if (i0.lt.1.or.mask(i0,j).eq.0) then     
            nbt_lc(k) = nbt_lc(k) + 1
            lbt_lc(nbt_lc(k),k) = ipt
         else
            nrt_lc(k) = nrt_lc(k) + 1
            lrt_lc(nrt_lc(k),k) = ipt
         endif
         it_l(i,j,k) = i2   ! i2 value of takeoff point
         at_l(i,j,k) = alf        ! distance from i2 to takeoff point
         dt_l(i,j,k)  = 1.
   10    continue
      enddo

c     find params for right characteristics of interior points
      nrt_rc(k) = 0
      nbt_rc(k) = 0
      not_rc(k) = 0
      do ix = 1, nptix
         ipt = lptix(ix)
         i   = ival(ipt)
         j   = jval(ipt)
         do is = 1 , idtx
            mm = 0
            if (i+is+1.le.nx) mm = mask(i+is+1,j)
            if (mm.eq.0) then  !characteristic intersects boundary
               not_rc(k) = not_rc(k) + 1
               lot_rc(not_rc(k),k) = ipt
               it_r(i,j,k) = i + is  ! i value of boundary point
               dt_r(i,j,k) = is*hx/dt
               goto 20
            endif
         enddo
c            characteristic did not intersect boundary
c               and takeoff point is in the interval:  [i1,i2)
         i1 = i + idtx
         i2 = i1 + 1
         i3 = i1 + 2
         if (i3.gt.nx.or.mask(i3,j).eq.0) then     
            nbt_rc(k) = nbt_rc(k) + 1
            lbt_rc(nbt_rc(k),k) = ipt
         else
            nrt_rc(k) = nrt_rc(k) + 1
            lrt_rc(nrt_rc(k),k) = ipt
         endif
         it_r(i,j,k) = i1
         at_r(i,j,k) = alf        ! distance from i1 to takeoff point
         dt_r(i,j,k)  = 1.
   20    continue
      enddo

c     find params for left characteristics of eastern boundary
      nrt_eb(k) = 0
      nbt_eb(k) = 0
      not_eb(k) = 0
      do ix = 1, npt_e
         ipt = lpt_e(ix)
         i   = ival(ipt)
         j   = jval(ipt)
         do is = 1 , idtx
            mm = 0
            if (i-is-1.ge.1) mm = mask(i-is-1,j)
            if (mm.eq.0) then  !characteristic intersects boundary
               not_eb(k) = not_eb(k) + 1
               lot_eb(not_eb(k),k) = ipt
               it_l(i,j,k) = i - is  ! i value of boundary point
               dt_l(i,j,k) = is*hx/dt
               goto 30
            endif
         enddo
c            characteristic did not intersect boundary
c               and takeoff point is in the interval:  (i1,i2]
         i2 = i - idtx
         i1 = i2 - 1
         i0 = i2 - 2
         i3 = i2 + 1
         if (i0.lt.1.or.i3.gt.nx.or.mask(i0,j)*mask(i3,j).eq.0) then     
            nbt_eb(k) = nbt_eb(k) + 1
            lbt_eb(nbt_eb(k),k) = ipt
         else
            nrt_eb(k) = nrt_eb(k) + 1
            lrt_eb(nrt_eb(k),k) = ipt
         endif
         it_l(i,j,k) = i2
         at_l(i,j,k) = alf        ! distance from i2 to takeoff point
         dt_l(i,j,k)  = 1.
   30    continue
      enddo


c     find params for right characteristic of western boundary points
      nrt_wb(k) = 0
      nbt_wb(k) = 0
      not_wb(k) = 0
      do ix = 1, npt_w
         ipt = lpt_w(ix)
         i   = ival(ipt)
         j   = jval(ipt)
         do is = 1 , idtx
            mm = 0
            if (i+is+1.le.nx) mm = mask(i+is+1,j)
            if (mm.eq.0) then  !characteristic intersects boundary
               not_wb(k) = not_wb(k) + 1
               lot_wb(not_wb(k),k) = ipt
c               print*,'ot_wb',not_wb(k),ipt
               it_r(i,j,k) = i + is  ! i value of boundary point
               dt_r(i,j,k) = is*hx/dt
               goto 40
            endif
         enddo
c            characteristic did not intersect boundary
c               and takeoff point is in the interval:  [i1,i2)
         i1 = i + idtx
         i2 = i1 + 1
         i3 = i1 + 2
         i0 = i1 - 1
         if (i0.lt.1.or.i3.gt.nx.or.mask(i0,j)*mask(i3,j).eq.0) then     
            nbt_wb(k) = nbt_wb(k) + 1
            lbt_wb(nbt_wb(k),k) = ipt
c               print*,'bt_wb',nbt_wb(k),ipt
         else
            nrt_wb(k) = nrt_wb(k) + 1
            lrt_wb(nrt_wb(k),k) = ipt
c               print*,'rt_wb',nrt_wb(k),ipt
         endif
         it_r(i,j,k) = i1
         at_r(i,j,k) = alf        ! distance from i1 to takeoff point
         dt_r(i,j,k)  = 1.
   40    continue
      enddo

c     make a list of all points within nshapd of the boundary
  200 continue

      if (k.eq.1) then

      do j = 1, ny
         do i = 1, nx
            weight_x(j,i) = 0.
            weight_y(j,i) = 0.
         enddo
      enddo

      n_yclose = 0
      do ind = 1, nptiy
         ipt = lptiy(ind)
         i = ival(ipt)
         j = jval(ipt)
         ndx = nx+ny
         ndy = nx+ny
         do in = 1, npt_e
            idpt = lpt_e(in)
            id = ival(idpt)
            jd = jval(idpt)
            if (ifper) then
               if (id.lt.nx) ndx = min(ndx,abs(id-i)+abs(jd-j))
            else
               ndx = min(ndx,abs(id-i)+abs(jd-j))
            endif
         enddo
         do in = 1, npt_w
            idpt = lpt_w(in)
            id = ival(idpt)
            jd = jval(idpt)
            if (ifper) then
               if (id.gt.1) ndx = min(ndx,abs(id-i)+abs(jd-j))
            else
               ndx = min(ndx,abs(id-i)+abs(jd-j))
            endif
         enddo
         do in = 1, npt_n
            idpt = lpt_n(in)
            id = ival(idpt)
            jd = jval(idpt)
            ndy = min(ndy,abs(id-i)+abs(jd-j))
         enddo
         do in = 1, npt_s
            idpt = lpt_s(in)
            id = ival(idpt)
            jd = jval(idpt)
            ndy = min(ndy,abs(id-i)+abs(jd-j))
         enddo
         if (ndx.le.nshapxd.or.ndy.le.nshapyd) then
             n_yclose = n_yclose+1
             l_yclose(n_yclose) = ipt
         endif
      enddo

      n_xclose = 0
      do ind = 1, nptix
         ipt = lptix(ind)
         i = ival(ipt)
         j = jval(ipt)
         ndx = nx+ny
         ndy = nx+ny
         do in = 1, npt_e
            idpt = lpt_e(in)
            id = ival(idpt)
            jd = jval(idpt)
            if (ifper) then
               if (id.lt.nx) ndx = min(ndx,abs(id-i)+abs(jd-j))
            else
               ndx = min(ndx,abs(id-i)+abs(jd-j))
            endif
         enddo
         do in = 1, npt_w
            idpt = lpt_w(in)
            id = ival(idpt)
            jd = jval(idpt)
            if (ifper) then
               if (id.gt.1) ndx = min(ndx,abs(id-i)+abs(jd-j))
            else
               ndx = min(ndx,abs(id-i)+abs(jd-j))
            endif
         enddo
         do in = 1, npt_n
            idpt = lpt_n(in)
            id = ival(idpt)
            jd = jval(idpt)
            ndy = min(ndy,abs(id-i)+abs(jd-j))
         enddo
         do in = 1, npt_s
            idpt = lpt_s(in)
            id = ival(idpt)
            jd = jval(idpt)
            ndy = min(ndy,abs(id-i)+abs(jd-j))
         enddo
         if (ndx.le.nshapxd.or.ndy.le.nshapyd) then
             n_xclose = n_xclose+1
             l_xclose(n_xclose) = ipt
         endif
      enddo
      do ind = 1, n_yclose
         ipt = l_yclose(ind)
         i = ival(ipt)
         j = jval(ipt)
         weight_y(j,i) = 1.
      enddo
      do ind = 1, n_xclose
         ipt = l_xclose(ind)
         i = ival(ipt)
         j = jval(ipt)
         weight_x(j,i) = 1.
      enddo

      endif

      return
      end


c----------------------------------------------------------------------|
      subroutine trisolv(np,rhs,v,a,b,c,bcleft,bcright)
      implicit none
      include 'om_para.h'
c----------------------------------------------------------------------|
c       tridiagonal solve for the system:
c   a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = rhs(i)
c       where  v(1) and v(n) are known
c       
c       on entry, rhs contains the righthand side
c       on exit , v contains the solution 
c----------------------------------------------------------------------|
      
      integer np,i,im
      real p(NYMAX),q(NYMAX),denom
      real a(np), b(np), c(np), v(np), rhs(np)
      real bcleft,bcright,pm,qm,am,vp
      
      v(1)  = bcleft
      vp = bcright
      
      pm = 0.0
      qm = v(1)
      
      do i = 2, np - 1
         im = i-1
         am = a(i)
         denom = 1./(am * pm + b(i))
         pm = -c(i) * denom
         qm = (rhs(i) - am * qm) * denom
         p(i) = pm
         q(i) = qm
      enddo
      
      do i = np - 1, 2, -1
         vp = p(i) * vp + q(i) 
         v(i) = vp
      enddo
      
      return
      end
      
c----------------------------------------------------------------------
      subroutine make_rhs4uh_cons(k,nxr,nyr,hm,um,vm,vmy,F,Q)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr
     
      integer i,j,in,ipt,i1,i2,ins,ind,ib
      real fc,alf,d1,db,ub,F(nyr,nxr),Q(nyr,nxr),d2
     *   ,hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr),vmy(nyr,nxr)
      real uphp,yvmvyp,fqp,umhm,yvpvym,fqm

      real uc(4,3),ud(3),ub_new,tk

      tk = tau(k)

      do in = 1, nrt_lc(k)
c---double characteristic case: left characteristic, cubic interpolation
         ipt = lrt_lc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 1,4
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) =    F(j,ind)   +  Q(j,ind)
         enddo
         call cub_interp(alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = d1
         dn(j,i) = d1
      enddo

      do in = 1, nbt_lc(k)
c---double characteristic case: left characteristic, linear interpolation
         ipt = lbt_lc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 2,3
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) = F(j,ind)   + Q(j,ind)
         enddo
         call lin_interp(alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = d1
         dn(j,i) = d1
      enddo

      do in = 1, not_lc(k)
c---double characteristic case: left characteristic intersects boundary
         ipt = lot_lc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         ib  = it_l(i,j,k)      
         db  = dt_l(i,j,k)      
         uphp   = um(j,ib)    + hm(j,ib)
         yvmvyp = fc*vm(j,ib) - vmy(j,ib)
         fqp    =  F(j,ib)    + Q(j,ib)
         d1 = 0.50*(uphp + tk*db*(fqp + yvmvyp))
         cn(j,i) = d1
         dn(j,i) = d1
      enddo


      do in = 1, nrt_rc(k)
c---double characteristic case: right characteristic, cubic interpolation
         ipt = lrt_rc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 1,4
            ind = i1 + ins - 2
            uc(ins,1) =    um(j,ind) -  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call cub_interp(1.-alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = cn(j,i) + d1
         dn(j,i) = dn(j,i) - d1
      enddo

      do in = 1, nbt_rc(k)
c---double characteristic case: right characteristic, linear interpolation
         ipt = lbt_rc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 2,3
            ind = i1 + ins - 2
            uc(ins,1) =    um(j,ind) -  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call lin_interp(1.-alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = cn(j,i) + d1
         dn(j,i) = dn(j,i) - d1
      enddo

      do in = 1, not_rc(k)
c---double characteristic case: right characteristic intersects boundary
         ipt = lot_rc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         ib  = it_r(i,j,k)      
         db = dt_r(i,j,k)      
         umhm   = um(j,ib)    - hm(j,ib)
         yvpvym = fc*vm(j,ib) + vmy(j,ib) 
         fqm    = F(j,ib)    - Q(j,ib)
         d1 = 0.50*(umhm + tk*db*(fqm + yvpvym))
         cn(j,i) = cn(j,i) + d1
         dn(j,i) = dn(j,i) - d1
      enddo

      do in = 1, nrt_eb(k)
c---single characteristic case: left characteristic, cubic interp
         ipt = lrt_eb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 1,4
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) =   F(j,ind) + Q(j,ind)
         enddo
         call cub_interp(alf,uc,ud)
         cn(j,i) = 0.
         d1 = hm(j,i) + tk*(fc*vm(j,i)-vmy(j,i)+F(j,i)+Q(j,i))
         d2 = ud(1) + tk*(ud(2) + ud(3))
         dn(j,i) = 0.5*(d1+d2)
      enddo

      do in = 1, nbt_eb(k)
c---single characteristic case: left characteristic, linear interp
         ipt = lbt_eb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 2,3
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) =   F(j,ind) + Q(j,ind)
         enddo
         call lin_interp(alf,uc,ud)
         cn(j,i) = 0.
         d1 = hm(j,i) + tk*(fc*vm(j,i)-vmy(j,i)+F(j,i)+Q(j,i))
         d2 = ud(1) + tk*(ud(2) + ud(3))
         dn(j,i) = 0.5*(d1+d2)
      enddo
      
      do in = 1, nrt_wb(k)
c---single characteristic case: right characteristic, cubic interp
         ipt = lrt_wb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)

         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 1,4
            ind = i1 + ins - 2
            uc(ins,1) =    hm(j,ind) -  um(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call cub_interp(1.-alf,uc,ud)
         cn(j,i) = 0.
         d1 = hm(j,i) - tk*(fc*vm(j,i)+vmy(j,i)+F(j,i)-Q(j,i))
         d2 = ud(1) - tk*(ud(2) + ud(3))
         dn(j,i) = 0.5*(d1+d2)
      enddo

      do in = 1, nbt_wb(k)
c---single characteristic case: right characteristic, linear interp
         ipt = lbt_wb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 2,3
            ind = i1 + ins - 2
            uc(ins,1) =    hm(j,ind) -  um(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call lin_interp(1.-alf,uc,ud)
         cn(j,i) = 0.
         d1 = hm(j,i) - tk*(fc*vm(j,i)+vmy(j,i)+F(j,i)-Q(j,i))
         d2 = ud(1) - tk*(ud(2) + ud(3))
         dn(j,i) = 0.5*(d1+d2)
      enddo
      
      return
      end

c----------------------------------------------------------------------
      subroutine make_rhs4uh(k,nxr,nyr,hm,um,vm,vmy,F,Q)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr
     
      integer i,j,in,ipt,i1,i2,ins,ind,ib
      real fc,alf,d1,db,ub,F(nyr,nxr),Q(nyr,nxr),d2
     *   ,hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr),vmy(nyr,nxr)
      real uphp,yvmvyp,fqp,umhm,yvpvym,fqm

      real uc(4,3),ud(3),ub_new,tk

      tk = tau(k)

      do in = 1, nrt_lc(k)
c---double characteristic case: left characteristic, cubic interpolation
         ipt = lrt_lc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 1,4
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) =    F(j,ind)   +  Q(j,ind)
         enddo
         call cub_interp(alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = d1
         dn(j,i) = d1
      enddo

      do in = 1, nbt_lc(k)
c---double characteristic case: left characteristic, linear interpolation
         ipt = lbt_lc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 2,3
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) = F(j,ind)   + Q(j,ind)
         enddo
         call lin_interp(alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = d1
         dn(j,i) = d1
      enddo

      do in = 1, not_lc(k)
c---double characteristic case: left characteristic intersects boundary
         ipt = lot_lc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         ib  = it_l(i,j,k)      
         db  = dt_l(i,j,k)      
         uphp   = um(j,ib)    + hm(j,ib)
         yvmvyp = fc*vm(j,ib) - vmy(j,ib)
         fqp    =  F(j,ib)    + Q(j,ib)
         d1 = 0.50*(uphp + tk*db*(fqp + yvmvyp))
         cn(j,i) = d1
         dn(j,i) = d1
      enddo


      do in = 1, nrt_rc(k)
c---double characteristic case: right characteristic, cubic interpolation
         ipt = lrt_rc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 1,4
            ind = i1 + ins - 2
            uc(ins,1) =    um(j,ind) -  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call cub_interp(1.-alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = cn(j,i) + d1
         dn(j,i) = dn(j,i) - d1
      enddo

      do in = 1, nbt_rc(k)
c---double characteristic case: right characteristic, linear interpolation
         ipt = lbt_rc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 2,3
            ind = i1 + ins - 2
            uc(ins,1) =    um(j,ind) -  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call lin_interp(1.-alf,uc,ud)

         d1 = 0.50*(ud(1) + tk*(ud(2) + ud(3)))
         cn(j,i) = cn(j,i) + d1
         dn(j,i) = dn(j,i) - d1
      enddo

      do in = 1, not_rc(k)
c---double characteristic case: right characteristic intersects boundary
         ipt = lot_rc(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         ib  = it_r(i,j,k)      
         db = dt_r(i,j,k)      
         umhm   = um(j,ib)    - hm(j,ib)
         yvpvym = fc*vm(j,ib) + vmy(j,ib) 
         fqm    = F(j,ib)    - Q(j,ib)
         d1 = 0.50*(umhm + tk*db*(fqm + yvpvym))
         cn(j,i) = cn(j,i) + d1
         dn(j,i) = dn(j,i) - d1
      enddo

      do in = 1, nrt_eb(k)
c---single characteristic case: left characteristic, cubic interp
         ipt = lrt_eb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 1,4
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) =   F(j,ind) + Q(j,ind)
         enddo
         call cub_interp(alf,uc,ud)
         cn(j,i) = 0.
         dn(j,i) = (ud(1) + tk*(ud(2) + ud(3)))
      enddo

      do in = 1, nbt_eb(k)
c---single characteristic case: left characteristic, linear interp
         ipt = lbt_eb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_l(i,j,k)      
         i2  = it_l(i,j,k)      
         do ins = 2,3
            ind = i2 + ins - 3
            uc(ins,1) =    um(j,ind) +  hm(j,ind)
            uc(ins,2) = fc*vm(j,ind) - vmy(j,ind)
            uc(ins,3) =   F(j,ind) + Q(j,ind)
         enddo
         call lin_interp(alf,uc,ud)
         cn(j,i) = 0.
         dn(j,i) = (ud(1) + tk*(ud(2) + ud(3)))
      enddo
      
      do in = 1, nrt_wb(k)
c---single characteristic case: right characteristic, cubic interp
         ipt = lrt_wb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)

         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 1,4
            ind = i1 + ins - 2
            uc(ins,1) =    hm(j,ind) -  um(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call cub_interp(1.-alf,uc,ud)
         cn(j,i) = 0.
         dn(j,i) = (ud(1) - tk*(ud(2) + ud(3)))
      enddo

      do in = 1, nbt_wb(k)
c---single characteristic case: right characteristic, linear interp
         ipt = lbt_wb(in,k)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         alf = at_r(i,j,k)      
         i1  = it_r(i,j,k)      
         do ins = 2,3
            ind = i1 + ins - 2
            uc(ins,1) =    hm(j,ind) -  um(j,ind)
            uc(ins,2) = fc*vm(j,ind) + vmy(j,ind)
            uc(ins,3) =   F(j,ind) - Q(j,ind)
         enddo
         call lin_interp(1.-alf,uc,ud)
         cn(j,i) = 0.
         dn(j,i) = ud(1) - tk*(ud(2) + ud(3))
      enddo
      
      return
      end

c----------------------------------------------------------------------
      subroutine add2rhs4uh_cons(k,nxr,nyr,F,Q)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr

      integer in,ipt,i,j
      real d1,d2,db,dp,tk,F(nyr,nxr),Q(nyr,nxr),tko2

      tk = tau(k)
      tko2 = 0.5*tk

      do in = 1, nptix
         ipt = lptix(in)
         i = ival(ipt)
         j = jval(ipt)
         d1 = dt_l(i,j,k)
         d2 = dt_r(i,j,k)
         db = tko2*(d1+d2)
         dp = tko2*(d1-d2)
         cn(j,i) = cn(j,i) + db*F(j,i) + dp*Q(j,i)
         dn(j,i) = dn(j,i) + db*Q(j,i) + dp*F(j,i)
      enddo
      do in = 1, npt_w
         ipt = lpt_w(in)
         i = ival(ipt)
         j = jval(ipt)
         dn(j,i) = dn(j,i) + tk*Q(j,i)
      enddo
      do in = 1, npt_e
         ipt = lpt_e(in)
         i = ival(ipt)
         j = jval(ipt)
         dn(j,i) = dn(j,i) + tk*Q(j,i)
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine add2rhs4uh(k,nxr,nyr,F,Q)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr

      integer in,ipt,i,j
      real d1,d2,db,dp,tk,F(nyr,nxr),Q(nyr,nxr)

      tk = tau(k)

      do in = 1, nptix
         ipt = lptix(in)
         i = ival(ipt)
         j = jval(ipt)
         d1 = dt_l(i,j,k)
         d2 = dt_r(i,j,k)
         db = tk*(d1+d2)/2.0
         dp = tk*(d1-d2)/2.0
         cn(j,i) = cn(j,i) + db*F(j,i) + dp*Q(j,i)
         dn(j,i) = dn(j,i) + db*Q(j,i) + dp*F(j,i)
      enddo
      do in = 1, npt_w
         ipt = lpt_w(in)
         i = ival(ipt)
         j = jval(ipt)
         dn(j,i) = dn(j,i) + tk*(Q(j,i) - F(j,i))
      enddo
      do in = 1, npt_e
         ipt = lpt_e(in)
         i = ival(ipt)
         j = jval(ipt)
         dn(j,i) = dn(j,i) + tk*(F(j,i) + Q(j,i))
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine make_rhs4v(k,nxr,nyr,hm,um,vm,G)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr

      integer in,ipt,i,j,jp,jm
      real d1,d2,db,fc,tyk,hk,tk,tik,G(nyr,nxr)
     *   ,hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr)

      tk = tau(k)
      tik = taui(k)
      hk = hysq(k)
      tyk = thyi(k)
      do in = 1, nptiy
         ipt = lptiy(in)
         i = ival(ipt)
         j = jval(ipt)
         d1 = dt_l(i,j,k)
         d2 = dt_r(i,j,k)
         db = tk*(d1+d2)/2.0
         jp = j+1
         jm = j-1
         fc = fcor(j,k)
         rhs(j,i) = (cn(j,i)+um(j,i) ) * fc +
     *        (dn(jp,i)+hm(jp,i)-dn(jm,i)-hm(jm,i))*tyk
     *        - G(j,i) - tik * vm(j,i)
         rhs(j,i) = hk*rhs(j,i)/db
      enddo

      return
      end
      
c----------------------------------------------------------------------
      subroutine add2rhs4v(k,nxr,nyr,G)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr

      integer in,ipt,i,j,jp,jm
      real d1,d2,db,fc,tk,hk,G(nyr,nxr)

      tk = tau(k)
      hk = hysq(k)
      do in = 1, nptiy
         ipt = lptiy(in)
         i = ival(ipt)
         j = jval(ipt)
         d1 = dt_l(i,j,k)
         d2 = dt_r(i,j,k)
         db = tk*(d1+d2)/2.0
         rhs(j,i) = rhs(j,i) - hk*G(j,i)/db
      enddo

      return
      end
      
c----------------------------------------------------------------------
      subroutine solve_v (nx, ny, mx, my, rr, vm , rgamm, rgamc, rgamp)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------

      integer i, nx, ny, mx, my
      real bcleft,bcright,rr(my,mx),vm(ny,nx)
     *     ,rgamm(my,mx),rgamc(my,mx),rgamp(my,mx)

      do i = 1,nx
         bcleft = 0.0
         bcright = 0.0
         call trisolv(ny,rr(1,i),vm(1,i),
     *        rgamm(1,i),rgamc(1,i),rgamp(1,i),bcleft,bcright)
      end do

      return
      end
      
c----------------------------------------------------------------------
      subroutine solve_uh_cons(k,nxr,nyr,hm,um,vm,vmy)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr

      integer in,ipt,i,j
      real fc,d1,d2,db,dp,tfv,tvy,tk
     *   ,hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr),vmy(nyr,nxr),tko2

      tk = tau(k)
      tko2 = 0.5*tk
      do in = 1, nptix
         ipt = lptix(in)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         d1 = dt_l(i,j,k)
         d2 = dt_r(i,j,k)
         db = tko2*(d1+d2)
         dp = tko2*(d1-d2)
         tfv = fc*vm(j,i)
         tvy = vmy(j,i)
         um(j,i) =  db*tfv - dp*tvy + cn(j,i)
         hm(j,i) =  dp*tfv - db*tvy + dn(j,i) 
      enddo
      do in = 1, npt_w
         ipt = lpt_w(in)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         um(j,i)= 0.0
         hm(j,i)= dn(j,i) - vmy(j,i)*tk 
      enddo
      do in = 1, npt_e
         ipt = lpt_e(in)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         um(j,i)= 0.0
         hm(j,i)= dn(j,i) - vmy(j,i)*tk
      enddo
      
      return
      end

c----------------------------------------------------------------------
      subroutine solve_uh(k,nxr,nyr,hm,um,vm,vmy)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr

      integer in,ipt,i,j
      real fc,d1,d2,db,dp,tfv,tvy,tk
     *   ,hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr),vmy(nyr,nxr)

      tk = tau(k)
      do in = 1, nptix
         ipt = lptix(in)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         d1 = dt_l(i,j,k)
         d2 = dt_r(i,j,k)
         db = tk*(d1+d2)/2.0
         dp = tk*(d1-d2)/2.0
         tfv = fc*vm(j,i)
         tvy = vmy(j,i)
         um(j,i) =  db*tfv - dp*tvy + cn(j,i)
         hm(j,i) =  dp*tfv - db*tvy + dn(j,i) 
      enddo
      do in = 1, npt_w
         ipt = lpt_w(in)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         um(j,i)= 0.0
         hm(j,i)= dn(j,i) - (fc*vm(j,i)+ vmy(j,i))*tk 
      enddo
      do in = 1, npt_e
         ipt = lpt_e(in)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)
         um(j,i)= 0.0
         hm(j,i)= dn(j,i) + (fc*vm(j,i)- vmy(j,i))*tk
      enddo
      
      return
      end

c----------------------------------------------------------------------
      subroutine make_vy(key,tyk,nxr,nyr,vm,vmy)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nxr,nyr,key
      real tyk,vm(nyr,nxr),vmy(nyr,nxr)

      integer in,ipt,i,j

      do in = 1, nptiy
         ipt = lptiy(in)
         i = ival(ipt)
         j = jval(ipt)
         vmy(j,i)=  (vm(j+1,i) - vm(j-1,i))* tyk
      enddo
      if (key.eq.0) then
         do in = 1, npt_n
            ipt = lpt_n(in)
            i = ival(ipt)
            j = jval(ipt)
            vmy(j,i)= -vm(j-1,i)*2.0*tyk
         enddo
         do in = 1, npt_s
            ipt = lpt_s(in)
            i = ival(ipt)
            j = jval(ipt)
            vmy(j,i)= vm(j+1,i)*2.0*tyk
         enddo
      else
         do in = 1, npt_n
            ipt = lpt_n(in)
            i = ival(ipt)
            j = jval(ipt)
            vmy(j,i)= -vm(j-1,i)*tyk
         enddo
         do in = 1, npt_s
            ipt = lpt_s(in)
            i = ival(ipt)
            j = jval(ipt)
            vmy(j,i)= vm(j+1,i)*tyk
         enddo
      endif
      
      return
      end

c----------------------------------------------------------------------
      subroutine bcset(nxr,nyr,um,vm)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
*  no-slip boundary conditions
      integer nxr,nyr
      real um(nyr,nxr),vm(nyr,nxr)

      integer in,ipt,i,j

      do in = 1, npt_ew
         ipt = lpt_ew(in)
         i = ival(ipt)
         j = jval(ipt)
         um(j,i) = 0.
      enddo
      do in = 1, npt_ns
         ipt = lpt_ns(in)
         i = ival(ipt)
         j = jval(ipt)
         vm(j,i) = 0.
      enddo
         
      return
      end

c----------------------------------------------------------------------
      subroutine bcset2(nxr,nyr,um,vm)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
*  no-slip boundary conditions
      integer nxr,nyr
      real um(nxr,nyr),vm(nxr,nyr)

      integer in,ipt,i,j

      do in = 1, npt_ew
         ipt = lpt_ew(in)
         i = ival(ipt)
         j = jval(ipt)
         um(i,j) = 0.
      enddo
      do in = 1, npt_ns
         ipt = lpt_ns(in)
         i = ival(ipt)
         j = jval(ipt)
         vm(i,j) = 0.
      enddo
         
      return
      end

c----------------------------------------------------------------------
      subroutine lin_interp(alf,uc,u)
      implicit none
c----------------------------------------------------------------------

      real alf, uc(4,3),u(3)

      integer i

      do i = 1, 3
         u(i) = alf*uc(2,i) + (1.0-alf)*uc(3,i)
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine cub_interp(alf,uc,u)
      implicit none
c----------------------------------------------------------------------

      real alf, uc(4,3),u(3)
      
      integer i
      real p

      do i = 1,3
         uc(4,i) = uc(4,i) - uc(3,i)
         uc(3,i) = uc(3,i) - uc(2,i)
         uc(2,i) = uc(2,i) - uc(1,i)
         uc(4,i) = (uc(4,i) - uc(3,i))/2.0
         uc(3,i) = (uc(3,i) - uc(2,i))/2.0
         uc(4,i) = (uc(4,i) - uc(3,i))/3.0
         
         p    = uc(3,i) - alf*uc(4,i)
         p    = uc(2,i) + (1.0- alf)*p
         u(i) = uc(1,i) + (2.0- alf)*p
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine shap(key,nx,ny,w,coefx,coefy)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nx, ny, key
      real w(ny,nx),coefx,coefy

      real x(NYMAX,NXMAX)
      integer in, ipt, i, j
      
c  filter y direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = w(j,i)
         end do
      end do

      do in = 1, nptiy
         ipt = lptiy(in)
         i   = ival(ipt)
         j   = jval(ipt)
         w(j,i) = x(j,i) + coefy*0.250*(x(j+1,i)+x(j-1,i)-2.0*x(j,i))
      end do
      if (key.eq.0) then  
         do in = 1, npt_n
            ipt = lpt_n(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefy*0.50*(x(j-1,i)-x(j,i))
         end do
         do in = 1, npt_s
            ipt = lpt_s(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefy*0.50*(x(j+1,i)-x(j,i))
         end do
      elseif (key.eq.1) then  
         do in = 1, npt_n
            ipt = lpt_n(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefy*0.250*(x(j-1,i)-x(j,i))
         end do
         do in = 1, npt_s
            ipt = lpt_s(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefy*0.250*(x(j+1,i)-x(j,i))
         end do
      endif

c  filter x direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = w(j,i)
         end do
      end do

      do in = 1, nptix
         ipt = lptix(in)
         i   = ival(ipt)
         j   = jval(ipt)
         w(j,i) = x(j,i) + coefx*0.250*(x(j,i+1)+x(j,i-1)-2.0*x(j,i))
      end do
      if (key.eq.0) then  
         do in = 1, npt_e
            ipt = lpt_e(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefx*0.50*(x(j,i-1)-x(j,i))
         end do
         do in = 1, npt_w
            ipt = lpt_w(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefx*0.50*(x(j,i+1)-x(j,i))
         end do
      elseif (key.eq.1) then  
         do in = 1, npt_e
            ipt = lpt_e(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefx*0.250*(x(j,i-1)-x(j,i))
         end do
         do in = 1, npt_w
            ipt = lpt_w(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coefx*0.250*(x(j,i+1)-x(j,i))
         end do
      endif

      return
      end

c----------------------------------------------------------------------
      subroutine shap_cons(nx,ny,h,shaph_coef)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
c
c     conservative filter for full box discretization near boundaries
c
c----------------------------------------------------------------------
      integer nx, ny, key
      real h(ny,nx),shaph_coef

      real x(NYMAX,NXMAX),coef,wc,wp,wm,hc,hyp,hym,hxp,hxm
      integer in, ipt, i, j
      
c  filter y direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = h(j,i)
         end do
      end do
  
      coef = 0.125 * shaph_coef

      do in = 1, nptiy
         ipt = lptiy(in)
         i   = ival(ipt)
         j   = jval(ipt)
         wc  = weight_y(j,i)
         wp  = weight_y(j+1,i)+wc
         wm  = weight_y(j-1,i)+wc
         hc  = x(j,i)
         hyp = x(j+1,i) - hc
         hym = hc - x(j-1,i)
         h(j,i) = x(j,i) + coef*(wp*hyp - wm*hym)
      end do
      do in = 1, npt_n
         ipt = lpt_n(in)
         i = ival(ipt)
         j = jval(ipt)
         wc = weight_y(j,i)
         wm = weight_y(j-1,i)+wc
         hc  = x(j,i)
         hym = hc - x(j-1,i)
         h(j,i) = x(j,i) - coef*wm*hym
      end do
      do in = 1, npt_s
         ipt = lpt_s(in)
         i = ival(ipt)
         j = jval(ipt)
         wc  = weight_y(j,i)
         wp  = weight_y(j+1,i)+wc
         hc  = x(j,i)
         hyp = x(j+1,i) - hc
         h(j,i) = x(j,i) + coef*wp*hyp
      end do

c  filter x direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = h(j,i)
         end do
      end do

      do in = 1, nptix
         ipt = lptix(in)
         i   = ival(ipt)
         j   = jval(ipt)
         wc  = weight_x(j,i)
         wp  = weight_x(j,i+1)+wc
         wm  = weight_x(j,i-1)+wc
         hc  = x(j,i)
         hxp = x(j,i+1) - hc
         hxm = hc - x(j,i-1)
         h(j,i) = x(j,i) + coef*(wp*hxp - wm*hxm)
      end do
      do in = 1, npt_e
         ipt = lpt_e(in)
         i = ival(ipt)
         j = jval(ipt)
         wc  = weight_x(j,i)
         wm  = weight_x(j,i-1)+wc
         hc  = x(j,i)
         hxm = hc - x(j,i-1)
         h(j,i) = x(j,i) - coef*wm*hxm
      end do
      do in = 1, npt_w
         ipt = lpt_w(in)
         i = ival(ipt)
         j = jval(ipt)
         wc  = weight_x(j,i)
         wp  = weight_x(j,i+1)+wc
         hc  = x(j,i)
         hxp = x(j,i+1) - hc
         h(j,i) = x(j,i) + coef*wp*hxp
      end do

      return
      end

c----------------------------------------------------------------------
      subroutine filt (nstep, nx, ny, h, u, v, vy, rdt)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nstep, nx, ny, k
      real h(ny,nx), u(ny,nx), v(ny,nx), rdt, vy(ny,nx)

      real shapi_coef, shapb_coef, shaph_coef, shapi_coefx, shapi_coefy
      integer nshapxd, nshapyd, nshapb, nshap, nshapt
      common /filter/ nshap,nshapb,nshapxd,nshapyd,shapi_coef,shapb_coef
     *               , shaph_coef, nshapt, shapi_coefx, shapi_coefy

      integer in
      
      if (nshapb.gt.0) then
         do in = 1, nshapt 
            if ( mod(nstep,nshapb).eq.0) then
               call shap_cons(nx,ny,h,shaph_coef)
               call shap_bound(1,nx,ny,u,shapb_coef)
               call shap_bound(2,nx,ny,v,shapb_coef)
c               call shap_bound_cons(1,nx,ny,u,shapb_coef)
c               call shap_bound_cons(2,nx,ny,v,shapb_coef)
            endif
         enddo
      endif
      
      if (nshap.gt.0) then
         if ( mod(nstep,nshap).eq.0) then
            call shap(0,nx,ny,h,shapi_coefx,shapi_coefy)
            call shap(2,nx,ny,u,shapi_coefx,shapi_coefy)
            call shap(2,nx,ny,v,shapi_coefx,shapi_coefy)
         endif
      endif

      call Rayleigh_friction2(nx,ny,u,v,vy,h,rdt)

      return
      end

c----------------------------------------------------------------------
      subroutine Rayleigh_friction (nx, ny, u, v, h, rdt)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nx, ny
      real u(ny,nx), v(ny,nx), h(ny,nx), rdt

      integer ipt, i, j, k
      real fexp,fexpb

      fexp = exp(-rdt)

      do ipt = 1, npt
         i = ival(ipt)
         j = jval(ipt)
         u(j,i)  = fexp*u(j,i)
         v(j,i)  = fexp*v(j,i)
         h(j,i)  = fexp*h(j,i)
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine shap_bound_cons(key,nx,ny,w,coef)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer key, nx, ny
      real w(ny,nx),coef
      integer i, j, in, ipt

      real x(NYMAX,NXMAX),co

      co = coef*0.25

c  filter y direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = w(j,i)
         end do
      end do
      
      do in = 1, n_yclose
         ipt = l_yclose(in)
         i = ival(ipt)
         j = jval(ipt)
         w(j,i) = x(j,i) + co*(x(j+1,i)+x(j-1,i)-2.0*x(j,i))
      end do
      if (key.eq.2) then
         do in = 1, npt_ns
            ipt = lpt_ns(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = 0.
         enddo
      elseif (key.ne.4) then
         do in = 1, npt_n
            ipt = lpt_n(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + co*(x(j-1,i)-x(j,i))
         end do
         do in = 1, npt_s
            ipt = lpt_s(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + co*(x(j+1,i)-x(j,i))
         end do
      endif
      
c  filter x direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = w(j,i)
         end do
      end do
      
      do in = 1, n_xclose
         ipt = l_xclose(in)
         i = ival(ipt)
         j = jval(ipt)
         w(j,i) = x(j,i) + co*(x(j,i+1)+x(j,i-1)-2.0*x(j,i))
      end do
      if (key.eq.1) then
         do in = 1, npt_ew
            ipt = lpt_ew(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = 0.
         enddo
      elseif (key.ne.4) then
         do in = 1, npt_e
            ipt = lpt_e(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + co*(x(j,i-1)-x(j,i))
         end do
         do in = 1, npt_w
            ipt = lpt_w(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + co*(x(j,i+1)-x(j,i))
         end do
      endif
         
      return
      end

c----------------------------------------------------------------------
      subroutine shap_bound(key,nx,ny,w,coef)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer key, nx, ny
      real w(ny,nx),coef
      integer i, j, in, ipt

      real x(NYMAX,NXMAX),co

c  filter y direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = w(j,i)
         end do
      end do

      do in = 1, n_yclose
         ipt = l_yclose(in)
         i = ival(ipt)
         j = jval(ipt)
         w(j,i) = x(j,i) + coef*0.250*(x(j+1,i)+x(j-1,i)-2.0*x(j,i))
      end do
      if (key.eq.2) then
         do in = 1, npt_ns
            ipt = lpt_ns(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = 0.
         enddo
      elseif (key.ne.4) then
         do in = 1, npt_n
            ipt = lpt_n(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coef*0.50*(x(j-1,i)-x(j,i))
         end do
         do in = 1, npt_s
            ipt = lpt_s(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coef*0.50*(x(j+1,i)-x(j,i))
         end do
      endif

c  filter x direction:
      do i = 1,nx
         do j = 1,ny
            x(j,i) = w(j,i)
         end do
      end do

      do in = 1, n_xclose
         ipt = l_xclose(in)
         i = ival(ipt)
         j = jval(ipt)
         w(j,i) = x(j,i) + coef*0.250*(x(j,i+1)+x(j,i-1)-2.0*x(j,i))
      end do
      if (key.eq.1) then
         do in = 1, npt_ew
            ipt = lpt_ew(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = 0.
         enddo
      elseif (key.ne.4) then
         do in = 1, npt_e
            ipt = lpt_e(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coef*0.50*(x(j,i-1)-x(j,i))
         end do
         do in = 1, npt_w
            ipt = lpt_w(in)
            i = ival(ipt)
            j = jval(ipt)
            w(j,i) = x(j,i) + coef*0.50*(x(j,i+1)-x(j,i))
         end do
      endif

      return
      end

c----------------------------------------------------------------------
      subroutine filt_v (nstep, nx, ny, v, rdt)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nstep, nx, ny, k
      real v(ny,nx), rdt

      real shapi_coef, shapb_coef, shaph_coef, shapi_coefx, shapi_coefy
      integer nshapxd, nshapyd, nshapb, nshap, nshapt
      common /filter/ nshap,nshapb,nshapxd,nshapyd,shapi_coef,shapb_coef
     *               , shaph_coef, nshapt, shapi_coefx, shapi_coefy

      integer in
      
      if (nshapb.gt.0) then
         do in = 1, nshapt 
            if ( mod(nstep,nshapb).eq.0) then
               call shap_bound(2,nx,ny,v,shapb_coef)
            endif
         enddo
      endif
      
      if (nshap.gt.0) then
         if ( mod(nstep,nshap).eq.0) then
            call shap(2,nx,ny,v,shapi_coefx,shapi_coefy)
         endif
      endif

      return
      end

c----------------------------------------------------------------------
      subroutine filt_uh (nstep, nx, ny, h, u, v, vy, rdt)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nstep, nx, ny, k
      real h(ny,nx), u(ny,nx), v(ny,nx), rdt, vy(ny,nx)

      real shapi_coef, shapb_coef, shaph_coef, shapi_coefx, shapi_coefy
      integer nshapxd, nshapyd, nshapb, nshap, nshapt
      common /filter/ nshap,nshapb,nshapxd,nshapyd,shapi_coef,shapb_coef
     *               , shaph_coef, nshapt, shapi_coefx, shapi_coefy

      integer in
      
      if (nshapb.gt.0) then
         do in = 1, nshapt 
            if ( mod(nstep,nshapb).eq.0) then
               call shap_cons(nx,ny,h,shaph_coef)
               call shap_bound(1,nx,ny,u,shapb_coef)
c               call shap_bound_cons(1,nx,ny,u,shapb_coef)
c               call shap_bound_cons(2,nx,ny,v,shapb_coef)
            endif
         enddo
      endif
      
      if (nshap.gt.0) then
         if ( mod(nstep,nshap).eq.0) then
            call shap(0,nx,ny,h,shapi_coefx,shapi_coefy)
            call shap(2,nx,ny,u,shapi_coefx,shapi_coefy)
         endif
      endif

      call Rayleigh_friction2(nx,ny,u,v,vy,h,rdt)

      return
      end

c----------------------------------------------------------------------
      subroutine Rayleigh_friction2 (nx, ny, u, v, vy, h, rdt)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer nx, ny
      real u(ny,nx), v(ny,nx), h(ny,nx), rdt, vy(ny,nx)

      integer ipt, i, j, k
      real fexp,fexpb

      fexp = exp(-rdt)

      do ipt = 1, npt
         i = ival(ipt)
         j = jval(ipt)
         u(j,i)  = fexp*u(j,i)
         v(j,i)  = fexp*v(j,i)
         vy(j,i) = fexp*vy(j,i)
         h(j,i)  = fexp*h(j,i)
      enddo

      return
      end

c----------------------------------------------------------------------
      subroutine add_diff(k,nxr,nyr,hm,um,vm,F,G,Q)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
      integer k,nxr,nyr
      real F(nyr,nxr),G(nyr,nxr),Q(nyr,nxr)
     *     ,hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr)

      integer ik, ipt, i, j
      real coef, hx2, hy2, coefx, coefy
      
      coef = vnu(k)
      coefx = vnu(k)
c#why?      coefy = vnu(k)/10.
      coefy = vnu(k)
      if (coef.eq.0.) return

      hx2 = 1./hxsq(k)
      hy2 = 1./hysq(k)
         
      do ik = 1, nptix
         ipt = lptix(ik)
         i = ival(ipt)
         j = jval(ipt)
         F(j,i) = F(j,i) + coefx*(um(j,i+1)+ um(j,i-1)- 2.*um(j,i))*hx2
         G(j,i) = G(j,i) + coefx*(vm(j,i+1)+ vm(j,i-1)- 2.*vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)+ hm(j,i-1)- 2.*hm(j,i))*hx2
      enddo
      do ik = 1, npt_e
         ipt = lpt_e(ik)
         i = ival(ipt)
         j = jval(ipt)
         F(j,i) = F(j,i) + coefx*(um(j,i-1)- um(j,i))*hx2
         G(j,i) = G(j,i) + coefx*(vm(j,i-1)- vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i-1)- hm(j,i))*hx2
      enddo
      do ik = 1, npt_w
         ipt = lpt_w(ik)
         i = ival(ipt)
         j = jval(ipt)
         F(j,i) = F(j,i) + coefx*(um(j,i+1)- um(j,i))*hx2
         G(j,i) = G(j,i) + coefx*(vm(j,i+1)- vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)- hm(j,i))*hx2
      enddo
      
      do ik = 1, nptiy
         ipt = lptiy(ik)
         i = ival(ipt)
         j = jval(ipt)
         F(j,i) = F(j,i) + coefy*(um(j+1,i)+ um(j-1,i)- 2.*um(j,i))*hy2
         G(j,i) = G(j,i) + coefy*(vm(j+1,i)+ vm(j-1,i)- 2.*vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j+1,i)+ hm(j-1,i)- 2.*hm(j,i))*hy2
      enddo
      do ik = 1, npt_n
         ipt = lpt_n(ik)
         i = ival(ipt)
         j = jval(ipt)
         F(j,i) = F(j,i) + coefy*(um(j-1,i)- um(j,i))*hy2
         G(j,i) = G(j,i) + coefy*(vm(j-1,i)- vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j-1,i)- hm(j,i))*hy2
      enddo
      do ik = 1, npt_s
         ipt = lpt_s(ik)
         i = ival(ipt)
         j = jval(ipt)
         F(j,i) = F(j,i) + coefy*(um(j+1,i)- um(j,i))*hy2
         G(j,i) = G(j,i) + coefy*(vm(j+1,i)- vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j+1,i)- hm(j,i))*hy2
      enddo
     
      return
      end
