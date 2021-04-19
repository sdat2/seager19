c----------------------------------------------------------------
      subroutine init_data_out (tfile, dfile)
c----------------------------------------------------------------
      use wrap
      implicit none
      include 'om_wrap.h'
      character*(*) tfile, dfile
      real zz(1000)
      common /tios_id/ idh, idcph, idlight, idn2, idtp
      integer idh, idcph, idlight, idn2, idtp

      real pin(MAXNZ)
      common /vertt/ pin

      real rnan, alpha, beta, gam
      integer imap, imap2, i, j, k, idvar_tios, in
      integer id_g1,id_g2,id_g3,id_g4, tios_idvar, id_g5,id_g6,id_g7,id_g8
      save zz

      alpha = 0.
      beta = 0.
      gam = 0.

      if (nz.gt.551) print*,'need to increase size of zp in init_data_out'

      call tios_init (tfile, dfile, alpha, beta, gam)

      allocate (zim(nz))

      do i = 1, 1000
         zz(i) = real(i)
      enddo
      do i = 1, nz-1
         zim(i) = 0.5*(zin(i)+zin(i+1))
      enddo

      call tios_grid (id_g1, ny, nx, nmodes, yv, xv, zz)
      call tios_grid (id_g2, nx, ny, 1, xv, yv, zz)
      call tios_grid (id_g3, nx, ny, nz, xv, yv, zin)
      call tios_grid (id_g5, ny, nx, 1, yv, xv, zz)
      call tios_grid (id_g6, nx, ny, nz-1, xv, yv, zim)
      call tios_grid (id_g8, nx, ny, nmodes, xv, yv, zz)

      call tios_map  (imap, nx*ny, nxyc, iox)

      call tios_var (p, 'PRES'//CHAR(0), id_g3, 0)
      call tios_var (u, 'U_VEL'//CHAR(0), id_g3, 0)
      call tios_var (v, 'V_VEL'//CHAR(0), id_g3, 0)

      call tios_var (t, 'TEMP'//CHAR(0), id_g6, 0)

      idcph = idvar_tios ('CPH'//CHAR(0), id_g8, 0)
      idlight  = idvar_tios ('LIGHT'//CHAR(0), id_g8, 0)

c      idn2 = idvar_tios ('N2\0', id_g3, 0)
c      idtp = idvar_tios ('TEMP_PROF\0', id_g3, 0)

      call tios_var (hmod,'HMODEL'//CHAR(0), id_g2, 0)
      call tios_var (htherm,'HTHERM'//CHAR(0), id_g2, 0)
      call tios_var (tdeep, 'TDEEP'//CHAR(0), id_g2, 0)
      call tios_var (u1, 'U1'//CHAR(0), id_g2, 0)
      call tios_var (v1, 'V1'//CHAR(0), id_g2, 0)
      call tios_var (u_ekman, 'U_EKMAN'//CHAR(0), id_g2, 0)
      call tios_var (v_ekman, 'V_EKMAN'//CHAR(0), id_g2, 0)
      call tios_var (hmixx,'HMIX'//CHAR(0), id_g2, 0)

      call tios_var (sst,   'SST'//CHAR(0)  , id_g2, 0)
      call tios_var (qnet,  'QNET'//CHAR(0) , id_g2, 0)
      call tios_var (qfp,  'QFP'//CHAR(0) , id_g2, 0)
      call tios_var (qdiag,  'QFLX'//CHAR(0) , id_g2, 0)
      call tios_var (w1,  'W1'//CHAR(0) , id_g2, 0)
      call tios_var (sst_terms,  'UDTDX'//CHAR(0)  , id_g2, 0)
      call tios_var (sst_terms(1*nxy+1),  'VDTDY'//CHAR(0)  , id_g2, 0)
      call tios_var (sst_terms(2*nxy+1),  'UP_FLUX'//CHAR(0), id_g2, 0)
      call tios_var (sst_terms(3*nxy+1),  'QPRIME'//CHAR(0), id_g2, 0)
      call tios_var (sst_terms(4*nxy+1),  'QFC'//CHAR(0), id_g2, 0)
      call tios_var (sst_terms(5*nxy+1),  'TDMAX'//CHAR(0), id_g2, 0)
      call tios_var (sst_terms(6*nxy+1),  'WUP'//CHAR(0), id_g2, 0)

      call tios_var (sh,  'SH'//CHAR(0)  , id_g2, 0)
      call tios_var (rlh, 'RLH'//CHAR(0) , id_g2, 0)
      call tios_var (qlw, 'QLW'//CHAR(0) , id_g2, 0)
      call tios_var (slrt,'SLR'//CHAR(0) , id_g2, 0)
      call tios_var (cldt,'CLD'//CHAR(0) , id_g2, 0)
      call tios_var (ut,  'UWND'//CHAR(0), id_g2, 0)
      call tios_var (vt,  'VWND'//CHAR(0), id_g2, 0)
      call tios_var (wt,  'WNSP'//CHAR(0), id_g2, 0)
      call tios_var (att, 'ATEM'//CHAR(0), id_g2, 0)
      call tios_var (aht, 'AHUM'//CHAR(0), id_g2, 0)

      call tios_var (hm, 'H_MODE'//CHAR(0), id_g1, 0)
      call tios_var (um, 'U_MODE'//CHAR(0), id_g1, 0)
      call tios_var (vm, 'V_MODE'//CHAR(0), id_g1, 0)

      call tios_var (F, 'TAUX'//CHAR(0), id_g5, 0)
      call tios_var (G, 'TAUY'//CHAR(0), id_g5, 0)

      call tios_read

      rnan = -987654.0

      do k = 1, nx*ny
         htherm(k) = rnan
         hmod(k) = rnan
         tdeep(k) = rnan
      enddo

      do k = 1, nxy
         qnet(k) = rnan
         qfp(k) = rnan
         qdiag(k) = rnan
         sh(k) = rnan
         rlh(k) = rnan
         qlw(k) = rnan
         slrt(k) = rnan
         cldt(k) = rnan
         ut(k) = rnan
         vt(k) = rnan
         wt(k) = rnan
      enddo

      return
      end

c-----------------------------------
      subroutine data_out (time)
      use wrap
      include 'om_wrap.h'
c-----------------------------------
      real time

      common /tios_id/ idh, idcph, idlight, idn2, idtp
      integer idh, tios_putidvar, idcph, idlight, idn2, idtp, tios_putvar
      external make_huv, flip_cph, flip_n2, flip_temp, flip_light

      if (n_debug.ne.0) print*,'beginning of data_out'

c      CALL TIOS_CNTRL (1, 1)   ! to dump every time step

      if (n_debug.ne.0) print*,time

      if (tios_putvar (p, time, make_huv).eq.0) then
         iret = tios_putvar (t, time, make_huv)
      else
         iret = tios_putvar (t, time, 0)
      endif

      iret = tios_putvar (tp, time, 0)

      iret = tios_putvar (hmod, time, 0)
      iret = tios_putvar (htherm, time, 0)

      iret =  tios_putvar (u_ekman, time, 0)
      iret =  tios_putvar (v_ekman, time, 0)

      iret = tios_putidvar(idcph, tp, time, flip_cph)
c      iret = tios_putidvar(idlight, tp, time, flip_light)
c      iret = tios_putvar(tp, time, flip_cph)
c      iret = tios_putvar(tp, time, flip_light)

c      iret = tios_putidvar(idn2, tp, time, flip_n2)
c      iret = tios_putidvar(idtp, tp, time, flip_temp)

      iret =  tios_putvar (hmixx, time, 0)

      iret =  tios_putvar (sst, time, 0)

      iret =  tios_putvar (qnet, time, 0)
      iret =  tios_putvar (qfp, time, 0)
      iret =  tios_putvar (qdiag, time, 0)
      iret =  tios_putvar (sst_terms, time, 0)
      iret =  tios_putvar (sst_terms(nxym+1), time, 0)
      iret =  tios_putvar (sst_terms(2*nxym+1), time, 0)
      iret =  tios_putvar (sst_terms(3*nxym+1), time, 0)
      iret =  tios_putvar (sst_terms(4*nxym+1), time, 0)

      iret =  tios_putvar (w1, time, 0)

      
      if (n_debug.ne.0) print*,'before tios_putvar'
      iret =  tios_putvar (hm, time, 0)
      if (n_debug.ne.0) print*,'after tios_putvar'

      iret =  tios_putvar (F, time, 0)
      iret =  tios_putvar (G, time, 0)

      if (n_debug.ne.0) print*,'before tios_save'
      call tios_save
      if (n_debug.ne.0) print*,'after tios_save'

      return
      end


c-----------------------------------
      subroutine make_huv
c-----------------------------------
      use wrap
      include 'om_wrap.h'
      logical non_stable

      call sum_modes(nx,ny,nxyc,nz,nmodes,hscl,uscl,dscl,hvm,vvm,hcm,
     *               hm,um,vm,p,u,v,nmixx,u1,v1,t,zin,use_ekm)
c      iflag = non_stable(iout, nx, ny, nxyc,nz, u, v)

      return
      end

c-----------------------------------
      subroutine sum_modes(nx,ny,nxyc,nz,nmodes,hscl,uscl,dscl,hvm,vvm,hcm,
     *                     hm,um,vm,p,u,v,nmixx,u1,v1,t,zin,use_ekm)
c-----------------------------------
      include 'om_core.h'
      dimension hm(ny,nx,nmodes),um(ny,nx,nmodes),vm(ny,nx,nmodes)
     *         ,p(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz)
     *         ,hvm(nz,nmodes,nx,ny),vvm(nz,nmodes,nx,ny),hcm(nz,nmodes,nx,ny)
     *         ,hscl(nmodes),uscl(nmodes),dscl(nmodes,nx,ny)
     *         ,u1(nx,ny),v1(nx,ny),t(nx,ny,nz),zin(nz),nmixx(nx,ny)
      logical use_ekm

      mode = 1 
      p =  -987654
      u =   -987654
      v =   -987654

      do k = 1, nz
         do n = 1, nxyc
            i = ival(n) 
            j = jval(n) 
            ah   = dscl(mode,i,j)*hvm(k,mode,i,j)
            h_ah = hscl(mode)*ah
            u_ah = uscl(mode)*ah
            p(i,j,k) = h_ah*hm(j,i,mode)
            u(i,j,k) = u_ah*um(j,i,mode)
            v(i,j,k) = u_ah*vm(j,i,mode)
         enddo
      enddo

      do mode = 2, nmodes
         do k = 1, nz
            do n = 1, nxyc
               i = ival(n) 
               j = jval(n) 
               ah   = dscl(mode,i,j)*hvm(k,mode,i,j)
               h_ah = hscl(mode)*ah
               u_ah = uscl(mode)*ah
               p(i,j,k) = p(i,j,k) + h_ah*hm(j,i,mode)
               u(i,j,k) = u(i,j,k) + u_ah*um(j,i,mode)
               v(i,j,k) = v(i,j,k) + u_ah*vm(j,i,mode)
            enddo
         enddo
      enddo

      if (use_ekm) then
      do n = 1, nxyc
         i = ival(n) 
         j = jval(n) 
         kmax = nmixx(i,j)
         do k = 1, kmax
            u(i,j,k) = u1(i,j)
            v(i,j,k) = v1(i,j)
         enddo
      enddo
      endif

      do k = 1, nz-1
         do n = 1, nxyc
            i = ival(n) 
            j = jval(n) 
            kp = k+1
            t(i,j,k) = -(p(i,j,kp) - p(i,j,k))/(zin(kp)-zin(k))
         enddo
      enddo

      return
      end

c     --------------------------------------------------------
      logical function non_stable (iout, nx, ny, nxyc ,nz, u, v)
c     --------------------------------------------------------
      include 'om_core.h'
c     check to see that velocities are not bizarre
      integer iout, nx, ny, nz, iox(1),nxyc
      real u(nx,ny,nz),v(nx,ny,nz)

      integer i, icheck
c
      non_stable = .false.
      icheck = 0
      do n = 1, nxyc
         i = ival(n) 
         j = jval(n) 
         if (u(i,j,1)**2+v(i,j,1)**2.gt.400.) icheck = icheck + 1
      enddo
      
      non_stable = (icheck .ne. 0)
      if (non_stable) then
         write (iout, *) 'Number of illegal points =', icheck 
         do n = 1, nxyc
            i = ival(n) 
            j = jval(n) 
            if (u(i,j,1)**2 + v(i,j,1)**2 .gt. 400.) then
               write (iout, 12) i, j
            endif
         enddo
      endif
   12 format ('i,j=',2i8,'  velocity is greater than 20 m/s')

      return
      end


c-----------------------------------
      subroutine flip_cph
c-----------------------------------
      use wrap
      include 'om_wrap.h'

      call transpose(nx,ny,nmodes,nmodes,cph,tp)

      return
      end

c-----------------------------------
      subroutine flip_light
c-----------------------------------
      use wrap
      include 'om_wrap.h'

      call transpose(nx,ny,nmodes,nmodes,cl,tp)

      return
      end

c-----------------------------------
      subroutine flip_n2
c-----------------------------------
      use wrap
      include 'om_wrap.h'

      call transpose(nx,ny,MAXNZ,nz,anbv2,tp)

      return
      end

c-----------------------------------
      subroutine flip_temp
c-----------------------------------
      use wrap
      include 'om_wrap.h'

      call transpose(nx,ny,MAXNZ,nz,temp,tp)

      return
      end

c-----------------------------------
      subroutine transpose(nx,ny,nmodes,nz,cph,tp)
c-----------------------------------
      real cph(nmodes,nx,ny),tp(nx,ny,nz)


      do mode = 1, nz
         do j = 1, ny
            do i = 1, nx
               tp(i,j,mode) = cph(mode,i,j)
            enddo
         enddo
      enddo

      return
      end


