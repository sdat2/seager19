c----------------------------------------------------
      subroutine vdens_init0(nx,ny,npt,mmz,nin,pin,anbv2,temp,hsfc,hmixx,
     *              z_cut,zcut) 
      use wdata
      implicit none
      integer nx, ny, mmz, nin, npt
      real pin(mmz), anbv2(mmz,nx*ny), temp(mmz,nx*ny), hsfc, hmixx(nx*ny)
      real zcut(nx*ny),z_cut
      include 'om_data.h'
c----------------------------------------------------

      integer i, j, ipt, k

      open (unit = 63, file = fbden(1:n_den))

      DO I = 1,mmz
         READ(63,*,END=200) PIN(I),ANBV2(I,1),TEMP(I,1)
      ENDDO

      print*,'trouble reading the density profile in ',fbden(1:n_den)
      STOP

  200 CONTINUE

      NIN = I - 1

      if (nin.gt.mmz) then
         print*,'must increase MAXNZ in om_para.h'
         stop
      endif

      hmixx(1) = hsfc
      zcut(1) = z_cut
      do i = 2, nx*ny
         zcut(i) = z_cut
         do k = 1, nin
            temp(k,i)  = temp(k,1)
            anbv2(k,i) = anbv2(k,1)
         enddo
         hmixx(i) = hsfc
      enddo

      return
      end

c----------------------------------------------------
      subroutine vdens_init(nx,ny,npt,mmz,nin,pin,temp,salt,
     *                      anbv2,hmixx,nmixx,tp,z_cut,zcut) 
      use wdata
      implicit none
      integer nx, ny, mmz, nin, nmixx(nx*ny), npt
      real pin(mmz), anbv2(mmz,nx*ny), temp(mmz,nx*ny), hmixx(nx*ny), tp(1)
      real salt(mmz,nx*ny),zcut(nx*ny),z_cut
      include 'om_data.h'
c----------------------------------------------------

      integer idf_d,ltem
      logical use_profile
      integer nvert_prof
      real xlon_prof, ylat_prof, hmix_param
      integer  ipt_prof
      common /vert_profile/ use_profile,xlon_prof, ylat_prof, nvert_prof,
     *                      ipt_prof, hmix_param

      real const, tmix,dzp,dzm,abv2p,abv2m,beta,talpha,grav,dtemp,smix
      integer ipt,k,kmix,i
      parameter (TALPHA = 2.55e-4)
      parameter (GRAV   = 9.780490)

      dtemp = hmix_param

      call odb_open(idf_t, fbpot(1:n_pot), 0)
      call odb_open(idf_s, fbsal(1:n_sal), 0)
      call data_on_model_grid(nx,ny,npt,idf_t,ltem,'potemp',tp,1)

      call odb_rddm(idf_t, 'Z', nzclm)
      allocate (zclm(nzclm))
      call odb_rdgr (idf_t, 'Z', nzclm, pin)
      nin = nzclm
      if (nin.gt.mmz) then
         print*,'must increase MAXNZ in om_para.h'
         stop
      endif

      call read_linz(nx,ny,idf_t,ltem,npt,mpack,nin,mmz,pin,temp,tp,'potemp')
      call read_linz(nx,ny,idf_s,ltem,npt,mpack,nin,mmz,pin,salt,tp,'salt')

      if (nvert_prof.eq.1) then
c        compute hmixx
         do ipt = 1, nx*ny
            pin(1) = 0.0   ! reset first depth to surface
            tmix = temp(1,ipt)
            smix = salt(1,ipt)
            zcut(ipt) = z_cut
            do k = 2, nin
               if (temp(k,ipt).lt.-100) then
                  temp(k,ipt) = temp(k,ipt) + 999.
                  zcut(ipt) = min(pin(k),zcut(ipt))
               endif
            enddo
            
c            write(17,*)
c            write(17,*)zcut(ipt)
c            write(17,*)
c            write(17,*)(temp(k,ipt),k=1,nin)


c            zcut(ipt) = max(z_cut,zcut(ipt))

            do k = 2, nin
               hmixx(ipt) = pin(k)
               kmix = k 
               tmix = tmix + temp(k,ipt) 
               smix = smix + salt(k,ipt) 
               if (temp(1,ipt).gt.dtemp+temp(k,ipt)) goto 100
            enddo
  100       continue

            nmixx(ipt) = kmix
            tmix = tmix + temp(kmix+1,ipt) 
            smix = smix + salt(kmix+1,ipt) 

c        reset temperatures in mixed layer to constant value
            tmix = tmix/(kmix+1)
            smix = smix/(kmix+1)
c            do k = 1, kmix+1
c               temp(k,ipt) = tmix
c               salt(k,ipt) = smix
c            enddo

c        compute N^2 = -g (drho/dz) / rho_0
            call BVFS(nin,pin,temp(1,ipt),salt(1,ipt),anbv2(1,ipt))
            if (mod(ipt,100).eq.1)then 
            write(16,*)ipt
            write(16,*)
     *         (pin(k),temp(k,ipt),salt(k,ipt),anbv2(k,ipt),k=1,nin)
            endif

         enddo

      elseif (nvert_prof.eq.2) then
   
         ipt = ipt_prof

c        compute hmixx
         pin(1) = 0.0           ! reset first depth to surface
         tmix = temp(1,ipt)
         smix = salt(1,ipt)
         do k = 2, nin
            hmixx(ipt) = pin(k)
            kmix = k 
            tmix = tmix + temp(k,ipt) 
            smix = smix + salt(k,ipt) 
            if (temp(1,ipt).gt.dtemp+temp(k,ipt)) goto 200
         enddo
  200    continue
         tmix = tmix + temp(kmix+1,ipt) 
         smix = smix + salt(kmix+1,ipt) 
         
c        reset temperatures in mixed layer to constant value
         tmix = tmix/(kmix+1)
         smix = smix/(kmix+1)
         do k = 1, kmix+1
            temp(k,ipt) = tmix
            salt(k,ipt) = smix
         enddo
         
c        compute N^2 = -g (drho/dz) / rho_0
         call BVFS(nin,pin,temp(1,ipt),salt(1,ipt),anbv2(1,ipt))
         
         do k = 1, nin
            print*,'reading vertical profile', PIN(k),temp(k,ipt),salt(k,ipt),ANBV2(k,ipt),TEMP(k,ipt)
         enddo

c        copy info to all other grid points
         do i = 1, ipt-1
            hmixx(i) = hmixx(ipt)
            do k = 1, nin
               temp(k,i) = temp(k,ipt)
               salt(k,i) = salt(k,ipt)
               anbv2(k,i)= anbv2(k,ipt)
            enddo
         enddo
         do i = ipt+1, nx*ny
            hmixx(i) = hmixx(ipt)
            do k = 1, nin
               temp(k,i) = temp(k,ipt)
               salt(k,i) = salt(k,ipt)
               anbv2(k,i)= anbv2(k,ipt)
            enddo
         enddo

      endif

      return
      end
      
c----------------------------------------------------
      subroutine forcing_init (nstart,nx,ny,xv,yv,area,F,G,Q,tp,hm,ym,scale_L)
      use wdata
      implicit none
      integer nstart, nx, ny
      real F(ny,nx),G(ny,nx),tp(1),xv(nx),yv(ny),Q(ny,nx),area,hm(ny,nx),
     *     ym(ny), scale_L  
      include 'om_data.h'    ! fbwnd,n_wnd,taux0,tauy0,taux1,tauy1
      include 'om_core.h'    ! npt,ival,jval
c----------------------------------------------------

c..........initalize the winds according to MTAU

      real TAUCON, API, TAUINV
      parameter (TAUCON = 1030., API = 3.14159265, TAUINV = 1./TAUCON)

      real bx,by,bx_width,by_width,mass_value
      common/ideal/bx,by,bx_width,by_width,mass_value
      common/winds/mtau,matau,tausc,atau,froude,itau_cos,drag_fac,mass
      integer mtau,matau,itau_cos,mass
      real tausc,atau,froude,drag_fac

      integer it1, it2, i, j, k, mcount, ic, ixy
      real taub, tauc, tscl, tau_conv, y1, y2, tmpx, q_per

      real enso_east, enso_west

      taub = drag_fac* tausc/TAUCON
      tauc = drag_fac* atau/TAUCON

      allocate (taux0(nx*ny))
      allocate (tauy0(nx*ny))
      allocate (taux1(nx*ny))
      allocate (tauy1(nx*ny))


      if     (mtau .eq. 1) then
c....................open the data files  (names from subroutine model_input)
         call odb_open(idf_tx, fbwnd(1:n_wnd)//'.x', 0)
         call odb_open(idf_ty, fbwnd(1:n_wnd)//'.y', 0)

c....................read in number of time slices
         call odb_rddm(idf_tx, 'T', ntau)

c....................allocate memory for ttau
         allocate (ttau(ntau))

c....................read in time grid
         call odb_rdgr(idf_tx, 'T', ntau, ttau)

c....................find indices of bracketing time slices, it1 and it2
         call it_catch (ntau, ttau, nstart, it1, it2, tscl)

c....................store it2 value for use in tau_lin
         itau = it2

c....................check whether on model or data grid
         call data_on_model_grid(nx,ny,npt,idf_tx, ltau, 'tau', tp, 1)

c....................read data for bracketing time slices, it1 and it2
         call read_zt (idf_tx, ltau, 1, it1, 'taux', tp, taux0) 
         call read_zt (idf_tx, ltau, 1, it2, 'taux', tp, taux1) 
         call read_zt (idf_ty, ltau, 1, it1, 'tauy', tp, tauy0) 
         call read_zt (idf_ty, ltau, 1, it2, 'tauy', tp, tauy1) 

         do ic = 1, npt
            i = ival(ic)
            j = jval(ic)
            ixy = (j-1)*nx+i
            taux0(ixy) = taub * taux0(ixy)
            taux1(ixy) = taub * taux1(ixy)
            tauy0(ixy) = tauc * tauy0(ixy)
            tauy1(ixy) = tauc * tauy1(ixy)
         enddo

         do ic = 1, npt
            i = ival(ic)
            j = jval(ic)
            ixy = (j-1)*nx+i
            F(j,i) = taux0(ixy) + tscl * (taux1(ixy) - taux0(ixy))
            G(j,i) = tauy0(ixy) + tscl * (tauy1(ixy) - tauy0(ixy))
         enddo

      elseif (mtau .eq. 3) then
c..........3 - annualy averaged climatology
         call odb_open(idf_tx, fbwnd(1:n_wnd)//'.x', 0)
         call odb_open(idf_ty, fbwnd(1:n_wnd)//'.y', 0)

         call odb_rddm(idf_tx, 'T', ntau)

         call data_on_model_grid(nx,ny,npt,idf_tx, ltau, 'tau', tp, 1)

         do k = 1, ntau
            call read_zt (idf_tx, ltau, 1, k, 'taux', tp, taux1) 
            call read_zt (idf_ty, ltau, 1, k, 'tauy', tp, tauy1) 

            do ic = 1, npt
               i = ival(ic)
               j = jval(ic)
               ixy = (j-1)*nx+i
               taux0(ixy) = taux0(ixy) + taux1(ixy)
               tauy0(ixy) = tauy0(ixy) + tauy1(ixy)
            enddo
         enddo

         do k = 1, npt
            i = ival(k)
            j = jval(k)
            ixy = (j-1)*nx+i
            F(j,i) = taub * taux0(ixy)/ntau
            G(j,i) = tauc * tauy0(ixy)/ntau
         enddo

      elseif (mtau .eq. 5) then
c..........5 - COSINE winds
         if     (itau_cos .eq. 0) then
            do j = 1, ny
               tp(j) = taub*cos(API*(yv(j))/80.)
            enddo
         elseif (itau_cos .eq. 1) then
            y1 = yv(1)
            y2 = yv(ny)
            do j = 1, ny
               tp(j) = taub*cos(2.*API*( (yv(j)-y1)/(y2-y1) - 0.5))
            enddo
         elseif (itau_cos .eq. 2) then
            y1 = yv(1)
            y2 = yv(ny)
            do j = 1, ny
               tp(j) = taub*cos(API*( (yv(j)-y1)/(y2-y1) - 0.5))
            enddo
         elseif (itau_cos .eq. 3) then
            do j = 1, ny
               tp(j) = taub
            enddo
         endif

         do k = 1, npt
            i = ival(k)
            j = jval(k)
            ixy = (j-1)*nx+i
            tmpx     = tp(j)
            F(j,i)   = tmpx
            taux0(ixy) = tmpx
            G(j,i)   = tauc
            tauy0(ixy) = tauc
         enddo
      elseif (mtau .eq. 6) then
c..........6 - ENSO winds

c  find the easternmost equatorial grid point
         enso_east = xv(1)
         enso_west = xv(nx)
         do k = 1, npt
            i = ival(k)
            j = jval(k)
            if (abs(yv(j)).lt.0.25.and.xv(i).gt.enso_east) then
               enso_east = xv(i)
               i_east = i
               j_east = j
            endif
            if (abs(yv(j)).lt.0.25.and.xv(i).lt.enso_west) then
               enso_west = xv(i)
            endif
         enddo

         enso_x1 = enso_west + enso_x1_fac*(enso_east-enso_west)
         enso_x2 = enso_west + enso_x2_fac*(enso_east-enso_west)

         print*,'enso forcing:', enso_west, enso_x1, enso_x2, enso_east
         print*, enso_mu, yv(ny/2), enso_k

         do j = 1, ny
c  must convert ym in meters to non-dimensional units!
            tp(j) = enso_k*exp(-enso_mu*(ym(j)/scale_L)**2)
         enddo

         do k = 1, npt
            i = ival(k)
            j = jval(k)
            ixy = (j-1)*nx+i
            if (xv(i).ge.enso_x1.and.xv(i).le.enso_x2) then
               j = jval(k)
               taux0(ixy) = tp(j)
            else
               taux0(ixy) = 0.
            endif
            F(j,i)   = 0.
            G(j,i)   = 0.
            tauy0(ixy) = 0.
         enddo

      endif

      if (mass .eq. 0) then
         do k = 1, npt
            i = ival(k)
            j = jval(k)
            Q(j,i) = 0.
         enddo
      elseif (mass .eq. 1) then
         mcount = 0
         do k = 1, npt
            i = ival(k)
            j = jval(k)
c            if (xv(i).gt.340.and.yv(j).gt.60) then
            if (xv(i).gt.260.and.yv(j).gt.50) then
c            if (xv(i).gt.260.and.yv(j).lt.00) then
c            if ((xv(i)-360)**2+(yv(j)+20)**2.lt.25) then
               mcount = mcount+1
            endif
         enddo
         q_per = mass_value*area/float(mcount)
         do k = 1, npt
            i = ival(k)
            j = jval(k)
               if (xv(i).gt.260.and.yv(j).gt.50) then
c               if (xv(i).gt.340.and.yv(j).gt.60) then
c            if ((xv(i)-375)**2+(yv(j)+20)**2.lt.25) then
                  Q(j,i) = q_per
c                  hm(j,i) = 1.e-4
               endif
         enddo
      elseif (mass .eq. 2) then
         mcount = 0
         do k = 1, npt
            i = ival(k)
            j = jval(k)
            if (abs(xv(i)-bx).lt.bx_width.and.abs(yv(j)-by).lt.by_width) then
               mcount = mcount+1
            endif
         enddo
         q_per = mass_value*area/float(mcount)
c         print*,bx,by,bx_width,by_width
c         print*,mcount,q_per,area
         do k = 1, npt
            i = ival(k)
            j = jval(k)
            if (abs(xv(i)-bx).lt.bx_width.and.abs(yv(j)-by).lt.by_width) then
                  Q(j,i) = q_per
c                  hm(j,i) = 1.e-4
               endif
         enddo
      endif

      return
      end

c     ------------------------------------------------------------------
      subroutine forcing (nstep,nx,ny,nxyc,xv,yv,F,G,Q,tp,hm)
      use wdata
      implicit none
      integer nstep,nx,ny,nxyc
      real F(ny,nx),G(ny,nx),tp(1),Q(ny,nx),xv(nx),yv(ny),hm(ny,nx)
      include 'om_data.h'  ! ntau,ttau,itau,idf_tx,ltau
      include 'om_core.h'  ! ival,jval,npt
c     ------------------------------------------------------------------

      common/winds/mtau,matau,tausc,atau,froude,itau_cos,drag_fac,mass
      integer mtau,matau,itau_cos,mass
      real tausc,atau,froude,drag_fac,str

      real TAUCON, API, TAUINV
      parameter (TAUCON = 1030., API = 3.14159265, TAUINV = 1./TAUCON)

      integer it1,it2,i,j,k,ixy
      real tscl, tau_conv, taub, tauc

      taub = drag_fac*tausc/TAUCON
      tauc = drag_fac*atau/TAUCON

      do k = 1, npt
         i = ival(k)
         j = jval(k)
         Q(j,i) = 0.
      enddo

      if (mtau.eq.6) then
c  make the wind depend on the value of h at the eastern equatorial grid point
         str = hm(j_east,i_east)
         str = tanh(enso_n*str)/max(enso_n,0.0001)
         do k = 1, npt
            i = ival(k)
            j = jval(k)
            ixy = (j-1)*nx+i
            F(j,i) = str*taux0(ixy)
         enddo

      else

      call it_catch (ntau, ttau, nstep, it1, it2, tscl)

      if (it2 .ne. itau) then
         itau = it2

         do i = 1, nx*ny
            taux0(i) = taux1(i)
            tauy0(i) = tauy1(i)
         enddo

         call read_zt (idf_tx, ltau, 1, it2, 'taux', tp, taux1) 
         call read_zt (idf_ty, ltau, 1, it2, 'tauy', tp, tauy1) 

         do k = 1, npt
            i = ival(k)
            j = jval(k)
            ixy = (j-1)*nx+i
            taux1(ixy) = taub * taux1(ixy)
            tauy1(ixy) = tauc * tauy1(ixy)
         enddo
      endif
      
      do k = 1, npt
         i = ival(k)
         j = jval(k)
         ixy = (j-1)*nx+i
         F(j,i) = taux0(ixy) + tscl * (taux1(ixy) - taux0(ixy))
         G(j,i) = tauy0(ixy) + tscl * (tauy1(ixy) - tauy0(ixy))
      enddo
      endif
      
      return
      end

c----------------------------------------------------------------------
      subroutine it_catch (NN, tt, nstep, it1, it2, tscl)
      implicit none
      include 'om_wrap.h' 
c----------------------------------------------------------------------
c.....Returns relative shift & indexes which are bracket nstep. 

c  INPUT:
c     NN
c     tt
c     nstep

c     enso_start
c     enso_scale

c  OUTPUT:
c     it1,it2
c     tscl
c-------------- 
      integer nn,nstep     
      real tt(nn)

      integer it1, it2
      real tscl
c----------------------------------------------------------------------
      real denso, tstep

      denso = enso_start + enso_scale * nstep

      tstep = tt(2) - tt(1)

      if (tt(nn)-tt(1)+tstep .eq. 12.) then
c.....Periodic Climatology Data
         denso = mod(denso, 12.)
         if (denso .lt. 0.) denso = denso + 12. 

         do it2 = 1, NN
            if (denso .lt. tt(it2)) goto 100 
         enddo
         denso = denso - 12.
         it2 = 1
  100    if (it2 .eq. 1) then 
            it1  = NN
            tscl = (12. - tt(NN) + denso)/(12. - tt(NN) + tt(1)) 
         else
            it1  = it2 - 1
            tscl = (denso - tt(it1))/(tt(it2) - tt(it1)) 
         endif
      else
c.....Non-periodic Data       
         do it2 = 1, NN
            if (denso .lt. tt(it2)) goto 200 
         enddo
         it2 = NN+1
         tscl = 0.
  200    if    (it2 .eq. 1) then
            it1  = 1
            tscl = 0.
         elseif(it2 .eq. NN+1) then
            it1  = NN
            it2  = NN
            tscl = 0.
         else
            it1  = it2 - 1
            tscl = (denso - tt(it1))/(tt(it2) - tt(it1)) 
         endif
      endif
      return
      end

c----------------------------------------------------------------------
      subroutine data_on_model_grid (nx,ny,npt,idf, lret, tag, tp, key)
      use wdata
      implicit none
      integer nx,ny,npt,idf, lret,key
      real tp(1)
      character*(*) tag
      include 'om_data.h'  
c----------------------------------------------------------------------

      logical grids_equiv

      lret = 1

      if (idatgr .eq. 0) then
c........check if data on the same grid as the model:
         if ( grids_equiv(tag,nx,ny,npt,idf, nx,ny,npt, 0,0, tp)) then
            lret = 0
c            write(6, *) tag, ' data on same grid as model '
         endif
         
c........check if data on the same grid as previous data:
      elseif ( grids_equiv(tag,nx,ny,npt,idf,mxp,myp,mpack,msx,msy,tp)) then 
         lret = 1
c         write(6, *) tag, ' data on same grid as previous data'
         
c........check if data on the same grid as the model:
      elseif ( grids_equiv(tag,nx,ny,npt,idf, nx,ny,npt, 0,0, tp)) then
         lret = 0
c         write(6, *) tag, ' data on model grid'
         
      else
         write(6, *) tag, 'Only one data GRID allowed! Stop.'
         stop
      endif
         
c      if (key.eq.0) idatgr = -1
      
      return
      end

c-------------------------------------------------------------------
      logical function grids_equiv (tag,nxp,nyp,npt,idf, 
     *                              kxp,kyp,kpack,ksx,ksy, tpp)
      use wdata
      use wrap
      implicit none
      integer nxp,nyp,npt,idf,kxp,kyp,kpack,ksx,ksy
      real tpp(nxp*nyp)
      character*(*) tag
      include 'om_wrap.h'
      include 'om_data.h'
c-------------------------------------------------------------------

      logical odb_ifatt, odb_ifvar
      integer m_p, m_x, m_y, m_sx, m_sy, m_seg
      real xer, yer

      m_x = 0
      call odb_rdgr(idf, 'X', m_x, tpp)
      xer = 1.e-6*(xv(nxp) - xv(1))  
      if (xv(1) .lt. tpp(1)-xer .or. xv(nxp) .gt. tpp(m_x)+xer) then
       write (6, *) tag,'!!! X grid of DATA must cover the model region'
       write (6, *) xv(1),tpp(1)-xer,xv(nxp),tpp(m_x)+xer
       stop
      endif

      m_sx = 0
      if ( odb_ifatt(idf, 'X', 'stretched') )
     *     call odb_getiattr(idf, 'X', 'ns', m_sx) 

      m_y = 0
      yer = 1.e-6*(yv(nyp) - yv(1))  
      call odb_rdgr(idf, 'Y', m_y, tpp(m_x+1))
      if (yv(1).lt. tpp(m_x+1)-yer .or. yv(nyp) .gt. tpp(m_x+m_y)+yer) then
       write (6, *) tag,'!!! Y grid of DATA must cover the model region'
       stop
      endif
      m_sy = 0

      if ( odb_ifatt(idf, 'Y', 'stretched') )
     *     call odb_getiattr(idf, 'Y', 'ns', m_sy) 

      if (odb_ifvar(idf,'MASK')) then
         m_p   = 0
         call odb_rddm(idf, 'NPACK', m_p)
         m_seg = 0
         call odb_rddm (idf, 'NMASK', m_seg)
         allocate (itp(m_seg+nxp+nyp))
         call odb_rdvar(idf, 'MASK', itp)
      else
         m_p   = m_x*m_y
         m_seg = 2
         allocate (itp(m_seg+nxp+nyp))
         itp(1) = 1
         itp(2) = m_p
      endif

      if (
c  m_p.ne.kpack .or. 
     *    m_x.ne.kxp   .or. m_y.ne.kyp .or.
     *    m_sx.ne.ksx  .or. m_sy.ne.ksy ) then

         grids_equiv = .FALSE.
         if (idatgr .eq. 0) then
            mpack = m_p
            mxp   = m_x
            msx   = m_sx
            myp   = m_y
            msy   = m_sy
            mseg  = m_seg

            call datagrid_memory(nxp,nyp,npt,tpp)
            call blin_indx(nxp,nyp,itp)
            call blin_coef(nxp,nyp,itp(mseg+1))
         else
            print*,tag,' new interpolation coefficients needed'
         endif
      else
         grids_equiv = .TRUE.
      endif

      deallocate (itp)

      return
      end

c-------------------------------------------------------
      function near(nx,ny,ixy,iseg)
      use wdata
      implicit none
      include 'om_data.h'
      integer ixy, iseg(mseg), nx, ny
c-------------------------------------------------------

      integer ix0, iy0, near, ix, iy, i, j
      real x0, y0, dmin, d

      iy0 = 1 + (ixy-1)/mxp
      ix0 = ixy - mxp*(iy0-1)
      x0 = xd(ix0)
      y0 = yd(iy0)

      near = iseg(1)
      iy = 1       + (near-1)/mxp
      ix = near - (iy-1)*mxp
      dmin = (xd(ix)-x0)**2 + (yd(iy)-y0)**2
      
      do j = 2, mseg
         i = iseg(j)
         iy = 1 + (i-1)/mxp
         ix = i - mxp*(iy-1)
         d = (xd(ix)-x0)**2 + (yd(iy)-y0)**2
         if (d .lt. dmin) then
            dmin = d
            near = i
         endif
      enddo

      return
      end

c--------------------------------------
      subroutine blin_coef(nxp,nyp,iw)
      use wrap
      use wdata
      implicit none
      integer nxp,nyp,iw(1)
      include 'om_wrap.h'
      include 'om_data.h'
c--------------------------------------

      integer i,j,k,i1,j1
      real fx,fy

      call bracket(mxp, xd, nxp, xv, iw)
      call bracket(myp, yd, nyp, yv, iw(nxp+1))

      do k = 1, nxyc
         j = 1 + (iox(k)-1)/nxp
         i = iox(k) - (j-1)*nxp
c.....find the i,j location for the four surrounding DATA grid points.
         i1 = iw(i)
         j1 = iw(nxp+j)
c.....x-sort index of four DATA points surrounding MODEL point (i,j).
         im2d(k) = i1 + (j1-1)*mxp
c.....find the interpolation ratios.
         fx = (xv(i)-xd(i1))/(xd(i1+1)-xd(i1))
         fy = (yv(j)-yd(j1))/(yd(j1+1)-yd(j1))


         blcf(k)        = (1.-fx)*(1.-fy)   
         blcf(k+nxyc)   = fx*(1.-fy)        
         blcf(k+2*nxyc) = (1.-fx)*fy        
         blcf(k+3*nxyc) = fx*fy             

      enddo

      return
      end

c     -----------------------------------------------------------------
      subroutine bracket(nx1,x1,nx2,x2,it)
c     -----------------------------------------------------------------
c     find the elements of x1 which bracket each element of x2.
c     returns it(i), for i=1,nx2 such that:
c          x1(it(i)) .le. x2(i)  .and.  x2(i) .le. x1(it(i)+1) 
c     nx1  = (input) length of x1.
c     x1   = (input) must have x1(i+1) .gt. x1(i), i=1,nx1-1.
c     nx2  = (input) length of x2.
c     x2   = (input) must have x2(i+1) .gt. x2(i), i=1,nx2-1.
c     it   = (output) nx2 indices of the lower side of the pair of
c            consecutive elements of x1 which bracket x2(i).
c
c     must input x1(1).le.x2(1) .and. x1(nx1).ge.x2(nx2).
c
      dimension x1(1),x2(1),it(1)
      i1 = 1
      do 20 i2=1,nx2
   10 if(x2(i2).ge.x1(i1) .and. x2(i2).le.x1(i1+1)) goto 20
      i1 = i1 + 1
      if(i1.lt.nx1-1) goto 10
   20 it(i2) = i1
      return
c     end of bracket.
      end

c---------------------------------      
      subroutine afill(n, a, v)
c---------------------------------      
      dimension a(1)
      do i = 1, n
         a(i) = v
      enddo
      return
      end

c------------------------------------------------------------------------
      subroutine blin_intr0(npt, mxp, ixd, im2d, blcf, fd, f)
      implicit none
c------------------------------------------------------------------------
      integer npt, mxp
      integer ixd(1), im2d(1)
      real blcf(npt,4), f(1), fd(1)

      integer i, i1

c  if uncompressed, then ixd(m)=m

      do i = 1, npt
         i1 = im2d(i)
         f(i) = blcf(i,1)*fd(ixd(i1))     + blcf(i,2)*fd(ixd(i1+1))
     *        + blcf(i,3)*fd(ixd(i1+mxp)) + blcf(i,4)*fd(ixd(i1+mxp+1))
      enddo

      return
      end

c------------------------------------------------------------------------
      subroutine blin_intr(npt, mxp, iox, ixd, im2d, blcf, fd, f)
      implicit none
c------------------------------------------------------------------------
      integer npt, mxp
      integer ixd(1), im2d(1), iox(1)
      real blcf(npt,4), f(1), fd(1)

      integer i, i1

c  if uncompressed, then ixd(m)=m

      do i = 1, npt
         i1 = im2d(i)
         f(iox(i)) = blcf(i,1)*fd(ixd(i1))     + blcf(i,2)*fd(ixd(i1+1))
     *        + blcf(i,3)*fd(ixd(i1+mxp)) + blcf(i,4)*fd(ixd(i1+mxp+1))
      enddo

      return
      end

c-------------------------------------------------------
      subroutine blin_indx (nx,ny,iseg)
      use wdata
      include 'om_data.h'
      integer iseg(2,mseg/2),nx,ny
c-------------------------------------------------------

      integer i,j,k,near,kn
      
c.....fill in ixd() according with the data compression iseg(2,mseg/2):      

      k = 0
      do i = 1, mseg/2
         do j = iseg(1,i), iseg(2,i)
            k = k + 1
            ixd(j) = k
         enddo
      enddo
c.....continue data to all points
      do k = 1, mxp*myp
         kn = near(nx,ny,k,iseg)
         if (ixd(k) .eq. 0) ixd(k) = ixd(kn)
      enddo

      return
      end


c-----------------------------------------------------------------
      subroutine read_zt (idf, key, iz, it, tag, ftmp, fdata) 
      use wdata
      use wrap
      implicit none
      include 'om_data.h'
      include 'om_wrap.h'
      integer idf, key, iz, it, i, ixy
      real ftmp(1), fdata(1)
      character*(*) tag
c-----------------------------------------------------------------
      logical odb_ifvar

      if (key .eq. 0) then      !! data on MODEL grid
         if (odb_ifvar(idf,'MASK')) then   !! compressed   (npt) grid
            call odb_rd1v3 (idf, iz, it, tag, fdata)
         else                              !! uncompressed (nx,ny) grid
            call odb_rd2v4 (idf, iz, it, tag, fdata)
         endif
      else                      !! data on a different grid
         if (odb_ifvar(idf,'MASK')) then   !! compressed   (npt) grid
            call odb_rd1v3 (idf, iz, it, tag, ftmp)
         else                              !! uncompressed (nx,ny) grid
            call odb_rd2v4 (idf, iz, it, tag, ftmp)
         endif
         call blin_intr(nxyc, mxp, iox, ixd, im2d, blcf, ftmp, fdata)
      endif


      return
      end

c-----------------------------------------------------------------
      subroutine read_xy (idf, it, tag, ftmp, fdata) 
      use wdata
      use wrap
      implicit none
      include 'om_data.h'
      include 'om_wrap.h'
      integer idf, it, i, ixy
      real ftmp(1), fdata(1)
      character*(*) tag
c-----------------------------------------------------------------
      logical odb_ifvar

      call odb_rd1v2 (idf, it, tag, fdata)

      return
      end

c----------------------------------------------------------------------------
      subroutine read_linz(nxp,nyp,idf,key,NPT,MPT,nzp,mmz,hdat,fdat,fvert,tag)
c----------------------------------------------------------------------------
      use wdata
      use wrap
      implicit none
      integer nxp,nyp,idf,key,npt,mpt,nzp,mmz
      integer i,k
      real hdat(npt,1),fdat(mmz,1), zvert(1),fvert(1)
      character*(*) tag
      logical odb_ifvar
      include 'om_wrap.h'  
      include 'om_data.h'  

      real, dimension(:,:), allocatable ::aa,bb

      if (key .eq. 0) then
         allocate (aa(npt,nzp))

         do k = 1, nzp
            if (odb_ifvar(idf,'MASK')) then   !! compressed   (npt) grid
               call odb_rd1v3(idf, k, 1, tag, aa(1,k))
            else                              !! uncompressed (nx,ny) grid
               call odb_rd2v4 (idf, k, 1, tag, aa(1,k))
            endif
         enddo
         
         do i = 1, nxp*nyp
            do k = 1, nzp
               fdat(k,i) = aa(1,k)
            enddo
         enddo

         do i = 1, npt
            do k = 1, nzp
               fdat(k,iox(i)) = aa(i,k)
            enddo
         enddo

         deallocate (aa)

      else
         allocate (bb(mpt,nzp))

         do k = 1, nzp
            if (odb_ifvar(idf,'MASK')) then
               call odb_rd1v3(idf, k, 1, tag, bb(1,k))
            else
               call odb_rd2v4 (idf, k, 1, tag, bb(1,k))
            endif

            call blin_intr0(npt, mxp, ixd, im2d, blcf, bb(1,k), fvert)

            do i = 1, nxp*nyp
               fdat(k,i) = bb(1,k)
            enddo
            do i = 1, npt
               fdat(k,iox(i)) = fvert(i)
            enddo

         enddo

         deallocate (bb)
      endif

      return
      end

c------------------------------------------------------------------
      subroutine BVFS(nin,zin,t,s,bf)
c------------------------------------------------------------------
c BVFS  Brunt-Vaisala Frequency Squared (N^2)
c===========================================================================
c
c DESCRIPTION:
c    Calculates Brunt-Vaisala Frequency squared (N^2) at the mid depths
c    from the equation,
c
c               -g      d(pdens)
c         N2 =  ----- x --------
c               pdens     d(z)
c
c INPUT:  (all must have same dimensions MxN)
c   S   = salinity    [psu]
c   T   = potential temperature 
c   z   = depth
c
c OUTPUT:
c   bf  = Brunt-Vaisala Frequency squared [s^-2]
c
c REFERENCES:
c   A.E. Gill 1982. p.54  eqn 3.7.15
c   "Atmosphere-Ocean Dynamics"
c   Academic Press: New York.  ISBN: 0-12-283522-0
c
c   Jackett, D.R. and McDougall, T.J. 1994.
c   Minimal adjustment of hydrographic properties to achieve static
c   stability.  submitted J.Atmos.Ocean.Tech.
c
c=========================================================================
      include 'om_para.h'
      
      integer nin
      real zin(nin),t(nin),s(nin),bf(nin)
      real*8 rh, rh0
      real dens(MAXNZ)

      grav = -9.8

      do k = 1, nin
          call dens_eos (0., t(k), s(k), rh0, rh)
          dens(k) = 1.e3 * real(rh - 1.d0)
      enddo

      dzp   = zin(2) - zin(1)
      abv2p = grav*(dens(1)-dens(2))/dzp
      do k = 2, nin - 1
         dzm   = dzp
         abv2m = abv2p
         dzp   = zin(k+1) - zin(k)
         abv2p = grav*(dens(k)-dens(k+1))/dzp
         beta  = dzp/(dzp+dzm)
         bf(k) = max(0.,(beta*abv2m + (1.-beta)*abv2p)/(1000.+dens(k)))
      enddo
      bf(1)   = 0.
      bf(nin) = 0.

      return
      end
            

c****************************************************************************
      subroutine dens_eos(pr, t, s, r0, rr)
c****************************************************************************
c  sub to compute density
c  calls sub 'sbulk', for secant bulk modulas
c
c  r0 is density at p = 0  -  returned in gr cm**3
c  rr is in situ density   -  returned
c

      implicit double precision (a-z)

      real*4  t, s, pr

c      dimension a(0:5), b(0:4), c(0:2)


       parameter
     1 (a0=999.842594d+00,a1=6.793952d-02,a2=-9.095290d-03,
     2   a3=1.001685d-04,a4=-1.120083d-06,a5=6.536332d-09,

     3  b0=8.24493d-01,b1=-4.0899d-03,b2=7.6438d-05,
     4  b3=-8.2467d-07,b4=5.3875d-09,

     5  c0=-5.72466d-03,c1=1.0227d-04,c2=-1.6546d-06,

     6  d=4.8314d-04)



      if (t.lt.-4.0 .or. t.gt.40.0) then
        r0 = -99.9
        rr = -99.9
        return
      else if (s.lt.0.0 .or. s.gt.42.0) then
        r0 = -99.9
        rr = -99.9
        return
      else if (pr.lt.0.0 .or. pr.gt.10000.0) then
        r0 = -99.9
        rr = -99.9
        return
      end if


      call sbulk(pr, t, s, kk)
c  secant bulk modulas (k) of seawater
c
c  density of snow
c
      rw = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t +a0
c
c  density at p = 0
c
      r0 = rw + s*((((b4*t + b3)*t + b2)*t + b1)*t + b0)
     *     + s*sqrt(s)*((c2*t + c1)*t + c0) + s*s*d
c
c  in situ density
c
      p = pr / 10.0
c  p is in bars
      rr = r0 / (1.d0 - p / kk)

      rr = rr / 1.d3
c  densities are returned in
      r0 = r0 / 1.d3
c  grams / cubic centimeter

      return
      end


c****************************************************************************
      subroutine sbulk(pr, t, s, kk)
c****************************************************************************
c
c  subroutines to calculate density, spec vol, secant bulk
c  modulas and alpha & beta
c  based on unesco 1981 report
c  equation of state for seawater - millero
c  programmer - c. greengrove, jan 1982
c  modified for hp - p mele, sep '82
c
c  range:
c    s =  0 to 42 (practical salinity)
c    t = -4 to 40 (c)
c    pr =  0 to 10000 (decibars)
c
c  other units:
c    density = kg/m3 **3
c    bulk deni mod.(k) = bars
c
c
c  kk is secant bulk modulas - returned
c

      implicit double precision (a-z)

      real*4  t, s, pr, s12
c  single precision


       parameter
     1 (e0=19652.21d+00,e1=148.4206d+00,e2=-2.327105d+00,
     2  e3=1.360477d-02,e4=-5.155288d-05,

     3  f0=54.6746d+00,f1=-.603459d+00,f2=1.09987d-02,f3=-6.167d-05,

     4  g0=7.944d-02,g1=1.6483d-02,g2=-5.3009d-04,

     5  h0=3.239908d+00,h1=1.43713d-03,h2=1.16092d-04,h3=-5.77905d-07,

     6  i0=2.2838d-03,i1=-1.0981d-05,i2=-1.6078d-06,

     7  j=1.91075d-04,

     8  k0=8.50935d-05,k1=-6.12293d-06,k2=5.2787d-08,

     9  m0=-9.9348d-07,m1=2.0816d-08,m2=9.1697d-10)



      if (t.lt.-4.0 .or. t.gt.40.0) then

c  range specifications
        kk = -99.9
        return
      else if (s.lt.0.0 .or. s.gt.42.0) then
        kk = -99.9
        return
      else if (pr.lt.0.0 .or. pr.gt.10000.0) then
        kk = -99.9
        return
      end if

      p = pr / 10.0
c  convert to bars

c  define sqrt(s)
      s12=sqrt(s)
c
c  secant bulk modulas (k) of seawater
c
c  pure water terms of sbm are w series
c

      kw = (((e4*t + e3)*t + e2)*t + e1)*t + e0
      aw = ((h3*t + h2)*t + h1)*t + h0
      bw = (k2*t + k1)*t + k0
c
c  coeff for final equation
c
      aa = aw + s*((i2*t + i1)*t + i0 + j*s12)
      bb = bw + s*((m2*t + m1)*t + m0)
c
c  sbm at p = 0 first term in the final eq
c
      ko = kw + s*(((f3*t + f2)*t + f1)*t + f0)
     *     + s*s12*((g2*t + g1)*t + g0)
c
c  final eq sbm
c
      kk = (bb*p + aa)*p + ko


      return
      end

