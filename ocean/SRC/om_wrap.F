#define ITYPE_RST 0
#define MACHINE_WORD 4

!     idig.f / "uphi" - model/
!--------------------------------------
      function idig(xxxxx, n)
      integer xxxxx, n, res
      integer idig
!
      res = xxxxx
      if (n .gt. 1) res = xxxxx / 10**(n-1)
      idig = mod (res, 10)

      return
      end


!----------------------------------------------------------------------
      subroutine logfile
      use wrap
      implicit none
      include 'om_wrap.h'
!----------------------------------------------------------------------
      common /run/ nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast
      integer nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast,lnblnk
      common /run2/ iyear, iday, isec
      integer iyear, iday, isec

      write(iout, *) 'Input  base name: <', fbi(1:n_in), '>'
      write(iout, *) 'Output base name: <', fbo(1:n_out), '>'

      write(iout,2) nstart,nlaststart,nskip,nsteps,iyear,iday,isec,DELT
    2 format(/1x,'parameter values used in this run:'//
     + ' nstart= ', i15 ,    '  nlaststart= ', i15 /
     + '  nskip= ', i15 ,    '      nsteps= ', i15 /
     + '   year= ', i15 ,    '         day= ', i15 /
     + '    sec= ', i15 /
     + '   delt= ', e15.7)
!
!     print out grid parameters.
!
      write(iout,6) X_MIN, X_MAX, Y_MIN, Y_MAX
    6 format(/
     + '   alat= ', e15.7 ,  '        blat= ', e15.7 /
     + '   alon= ', e15.7 ,  '        blon= ', e15.7)

      end

      subroutine segm_from_iox (npt, nseg, iox, iseg, mseg)
!------------------------------------------------------
      dimension iox(npt), iseg(2,mseg)

      ista   = iox(1)
      inext  = ista + 1
      nseg   = 0
      icount = 0

      do i = 2, npt
         icurr = iox(i)
         icount = icount + 1

         if (inext .ne. icurr) then
            nseg = nseg + 1
            iseg(1,nseg) = ista
            iseg(2,nseg) = ista + icount - 1
            ista  = icurr
            icount = 0
         endif

         inext = icurr + 1
      enddo

      nseg = nseg + 1
      iseg(1,nseg) = ista
      iseg(2,nseg) = ista + icount

      return
      end

!----------------------------------------------------------------------
      subroutine model_input (tenso,inp_modes,phasev,
     *             time_mul,run_time,save_time,gamma,rs,f1prime)
      use wrap
      use wdata
      implicit none
      include 'om_wrap.h'
      include 'om_data.h'
!----------------------------------------------------------------------
      real tenso, phasev(1), time_mul,run_time,save_time,gamma,rs,f1prime

      real shapi_coef, shapb_coef, shaph_coef, shapi_coefx, shapi_coefy
      integer nshapxd, nshapyd, nshapb, nshap, nshapt
      common /filter/ nshap,nshapb,nshapxd,nshapyd,shapi_coef,shapb_coef
     *               , shaph_coef, nshapt, shapi_coefx, shapi_coefy

      common /run/ nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast
      integer nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast

      common /run2/ iyear, iday, isec
      integer iyear, iday, isec

      common/winds/mtau,matau,tausc,atau,froude,itau_cos,drag_fac,mass
      integer mtau,matau,itau_cos,nret,mass
      real tausc,atau,froude,drag_fac
      real flt(100)
      save         flt

      real bx,by,bx_width,by_width,mass_value
      common/ideal/bx,by,bx_width,by_width,mass_value

      integer i_prof, j_prof

      integer NSPTSMAX
      PARAMETER(NSPTSMAX=20)
      REAL*8 spknots(NSPTSMAX), spcoeffs(NSPTSMAX) 
      REAL knots(NSPTSMAX), coeffs(NSPTSMAX) 
      REAL B0,B1,B2,HCUT,TCUT
      INTEGER NKNOTS,NSPK
      COMMON /SSTPARAM/ spknots,spcoeffs,B0,B1,B2,HCUT,TCUT,nspk,nknots

      real xk(3)

      logical inp_def
      integer inp_str, inp_int, inp_modes, i, j, id, im, ns, inp_rarr,nzz
      integer lnblnk
      real inp_flt, inp_days, date2enso, dxdeg, dydeg, cspeed, ek_timescale
      real xmin_t, ymin_t

      call inp_file(finp)

      if (inp_int('Trace'//CHAR(0),1) .ne. 0) 
     *     call inp_trace(finp(1:lnblnk(finp))//'.tr')

      n_debug = inp_int('Debug'//CHAR(0),0)    

      n_dir= inp_str('Output_dir'//CHAR(0),'output'//CHAR(0), fbdir)

      n_tios= inp_str('Base_file'//CHAR(0),
     *     fbdir(1:n_dir)//'/'//finp(1:lnblnk(finp))//CHAR(0),fbt)
      n_out = inp_str('Save_file'//CHAR(0), fbt(1:n_tios)//'.save'//CHAR(0), fbo)
      n_in  = inp_str('Restart_file'//CHAR(0), fbo(1:n_out)//CHAR(0), fbi)
      n_den = inp_str('Density_file'//CHAR(0),  'DATA/rzk.pro'//CHAR(0), fbden)
      n_pot = inp_str('Temp_file'//CHAR(0),  'temp_data'//CHAR(0), fbpot)
      n_sal = inp_str('Salt_file'//CHAR(0),  'salt_data'//CHAR(0), fbsal)


      n_wnd = inp_str('Wind_file'//CHAR(0),  'wind_data'//CHAR(0), fbwnd)

      n_sst = inp_str('SST_file'//CHAR(0),   'sst_data'//CHAR(0), fbsst)
      n_qflx = inp_str('QFLX_file'//CHAR(0),   'qflx_data'//CHAR(0), fbqflx)
      n_dQdT = inp_str('dQdT_file'//CHAR(0),   'dQdT_data'//CHAR(0), fbdQdT)
      n_dQdf = inp_str('dQdf_file'//CHAR(0),   'dQdf_data'//CHAR(0), fbdQdf)
      n_cld = inp_str('Cloud_file'//CHAR(0), 'cld_data'//CHAR(0), fbcld)
      n_slr = inp_str('Solar_file'//CHAR(0), 'slr_data'//CHAR(0), fbslr)
      n_wsp = inp_str('Wndspd_file'//CHAR(0),'wndspd_data'//CHAR(0), fwsp)
      n_uwd = inp_str('Uwnd_file'//CHAR(0),'uwnd_data'//CHAR(0), fuwd)
      n_vwd = inp_str('Vwnd_file'//CHAR(0),'vwnd_data'//CHAR(0), fvwd)
      n_ah  = inp_str('Airhum_file'//CHAR(0),'ahum_data'//CHAR(0), fah)
      n_at  = inp_str('Airtem_file'//CHAR(0),'atem_data'//CHAR(0), fat)

      irest = inp_int('Restart'//CHAR(0), 0)        ! 0 - start; 1 - restart
      ischeme = inp_int('Scheme'//CHAR(0), 1) 
!                       1 - INC ; 2 - Leapfrog
      idiag = inp_int('Qflux_diag'//CHAR(0), 0) 

      mtau  =  inp_int('Wind_forc'//CHAR(0),1)
      if (mtau .eq. 5) itau_cos = inp_int('Wind_cos'//CHAR(0), 0)
      mass  =  inp_int('Mass_forc'//CHAR(0),0)
      bx =  inp_flt('Mass_bx'//CHAR(0),-20.)
      by =  inp_flt('Mass_by'//CHAR(0), 55.)
      bx_width =  inp_flt('Mass_bx_width'//CHAR(0), 10.)
      by_width =  inp_flt('Mass_by_width'//CHAR(0), 10.)
      mass_value =  inp_flt('Mass_value'//CHAR(0),1.5e7)

      if (mtau .eq. 6) then
         enso_k = inp_flt('ENSO_k'//CHAR(0), 1.)
         enso_n = inp_flt('ENSO_n'//CHAR(0), .0001)
         enso_mu = inp_flt('ENSO_mu'//CHAR(0), 0.1)
         enso_x1_fac = inp_flt('ENSO_x1_fac'//CHAR(0), 0.35)
         enso_x2_fac = inp_flt('ENSO_x2_fac'//CHAR(0), 0.65)
      endif

      flt(1) = 1.
      flt(2) = 1.
      nret =  inp_rarr ('Wind_tauxy'//CHAR(0), 2, flt, flt)
      tausc = flt(1)
      atau  = flt(2)
      drag_fac =  inp_flt('Wind_drag_fac'//CHAR(0),1.8)

      use_sst = .false.
      msst =  inp_int('SST_forc'//CHAR(0),0)
      if (msst.eq.1) use_sst = .true.

      f1prime =  inp_flt('f1prime'//CHAR(0),0.000)

      trans_coef_sst = 0.
      rlx_time_sst = 86400.*inp_days ('Rlx_time_sst'//CHAR(0), 0.)
      if (rlx_time_sst.ne.0)
     *                 trans_coef_sst = 1./rlx_time_sst    

      use_ekm = .false.
      mekm =  inp_int('Ekman_forc'//CHAR(0),0)
      if (mekm.eq.1) use_ekm = .true.

      if (use_sst) then
         ek_timescale = inp_days ('Ekman_r'//CHAR(0), 2.)
         rs =  1./ek_timescale/24./3600.
         gamma =  inp_flt('SST_gamma'//CHAR(0),0.75)
         hsfc =  inp_flt('Ekman_hsfc'//CHAR(0),50.)
      endif

      ifper = inp_def('Periodic'//CHAR(0))

      n_mask = inp_str('Mask_file'//CHAR(0),  'no'//CHAR(0), fbmask)

      if (n_mask .eq. 2) then
         print*,'reading Grid_mask for land/sea mask'
         nx = inp_int('NX'//CHAR(0),20)
         ny = inp_int('NY'//CHAR(0),20)
         X_MIN = inp_flt('West'//CHAR(0),0.)
         X_MAX = inp_flt('East'//CHAR(0),40.)
         Y_MIN = inp_flt('South'//CHAR(0),20.)
         Y_MAX = inp_flt('North'//CHAR(0),60.)
         allocate (xv(nx))
         allocate (yv(ny))
         allocate ( mask(nx*ny))
         mask = -1
         call read_mask ('Grid_mask'//CHAR(0), NX, NY, mask)
      else
         print*,'reading file for land/sea mask: ', fbmask(1:n_mask)
         call odb_open(idf_mask, fbmask(1:n_mask), 0)
         call odb_rddm(idf_mask, 'X', nx)
         call odb_rddm(idf_mask, 'Y', ny)
         allocate (xv(nx))
         allocate (yv(ny))
         allocate ( mask(nx*ny))
         allocate (rmask(nx*ny))
         call odb_rdgr(idf_mask, 'X', nx, xv)
         X_MIN = xv(1)
         X_MAX = xv(nx)
         call odb_rdgr(idf_mask, 'Y', ny, yv)
         Y_MIN = yv(1)
         Y_MAX = yv(ny)
         rmask = -1.0
         call read_xy (idf_mask, 0, 'mask', tp, rmask)
         mask = int(rmask)
      endif

      NXY  = NX * NY

!      do j = ny, 1, -1
!         write(9,199)(mask((j-1)*nx+i),i=1,nx)
!      enddo
!  199 format(360i1)

      NZ = inp_int('NZ'//CHAR(0),1)
      allocate (zin(nz))
      zin = 0.

      nzz = inp_rarr('Z_profile'//CHAR(0), nz, zin, zin)

      if (nz .gt. nzz) then
         print*,'Z_profile or H_profile: not enough terms...'
         stop
      endif

      nz_mixed = 0
      if (use_ekm) then
         do i = 1, nz
            if (zin(i).le.hsfc) nz_mixed = i
         enddo
      endif

      nmodes = inp_int('NUMMODE'//CHAR(0),1)
      inp_modes = inp_int('HEQUIVF'//CHAR(0),0)
      use_profile = .false.

      if (inp_modes.eq.0) then
         use_profile = .true.
!                   TYPE: 0, ascii file for a single profile
!                   TYPE: 1, netcdf file (LEVITUS) with varying profiles
!                   TYPE: 2, netcdf file (LEVITUS) with constant profiles-
!                           specify PROFILE_LAT,PROFILE_LON to pick (lat,lon)
          nvert_prof = inp_int('VERTICAL_PROFILE_TYPE'//CHAR(0),0)
          hmix_param = inp_flt('VERTICAL_HMIX_PARAM'//CHAR(0),0.5)
          ylat_prof = inp_flt('PROFILE_LAT'//CHAR(0),0.)
          xlon_prof = inp_flt('PROFILE_LON'//CHAR(0),180.)
      else
         nmodes = 1
         NZ = 1
         c_1mode  = inp_flt('Cspeed_1mode'//CHAR(0),3.)
         H_1mode  = inp_flt('Depth_1mode'//CHAR(0),150.)
         H1_1mode = inp_flt('HSFC_1mode'//CHAR(0),50.)
         print*,'c_1mode, H_1mode, H1_1mode',c_1mode, H_1mode, H1_1mode

      endif

      n_zis = inp_str('ZIS'//CHAR(0),  'depth'//CHAR(0), fb_zis)
      z_top = inp_flt('ZTOP'//CHAR(0),0.)
      z_bot = inp_flt('ZBOT'//CHAR(0),5500.)
      z_cut = inp_flt('ZCUT'//CHAR(0),6000.)

      DZm    = ( Z_BOT - Z_TOP ) / NZ
         
      A_rayl = inp_flt('A_Rayleigh'//CHAR(0),  1.e-8)
      A_lapl = inp_flt('A_Laplacian'//CHAR(0), 4.e+3)
      nshap  =  inp_int('Shap_inter_step'//CHAR(0),0)
      nshapb =  inp_int('Shap_bound_step'//CHAR(0),1)
      nshapt =  inp_int('Shap_bound_repeat'//CHAR(0),1)
      nshapxd =  inp_int('Shap_bound_xdist'//CHAR(0),1)
      nshapyd =  inp_int('Shap_bound_ydist'//CHAR(0),1)
      shapi_coef =  inp_flt('Shap_inter_coef'//CHAR(0),1.)
      shapi_coefx =  inp_flt('Shap_inter_coefx'//CHAR(0),shapi_coef)
      shapi_coefy =  inp_flt('Shap_inter_coefy'//CHAR(0),shapi_coef)
      shapb_coef =  inp_flt('Shap_bound_coef'//CHAR(0),1.)
      shaph_coef =  inp_flt('Shap_h_coef'//CHAR(0),shapb_coef)

      delt = inp_days ('Time_step'//CHAR(0), 1./24.) *24.*3600.
      time_mul = inp_flt ('Time_mul'//CHAR(0), 1.)

      if (ischeme.eq.0) then
         if (time_mul.ne.1.) then
            print*,'SCHEME 0 is only coded for dt=dx, (must use Time_mul 1.)'
            stop
         endif
      endif

      xv = 0.
      yv = 0.

      dxdeg  = (X_MAX - X_MIN)/(nx - 1)   !  x spacing in degrees
      dydeg  = (Y_MAX - Y_MIN)/(ny - 1)   !  y spacing in degrees
 
      do i=1,nx
         xv(i)  = X_MIN + (i - 1) * dxdeg
      end do
      do j=1,ny
         yv(j)  = Y_MIN + (j - 1) * dydeg
      end do

      i_prof = 0
      do i=1,nx
         if (xv(i).lt.xlon_prof) i_prof = i
      enddo

      j_prof = 0
      do j=1,ny
         if (yv(j).lt.ylat_prof) j_prof = j
      enddo

      nstart = 1  

      run_time = inp_days('Run_time'//CHAR(0),   365.)  *24.*3600.
      save_time = inp_days('Save_step'//CHAR(0), 365.)  *24.*3600. 

      iyear = inp_int('Year'//CHAR(0),   2000)
      iday  = inp_int('Day'//CHAR(0),   0)
      isec  = inp_int('Sec'//CHAR(0),   0)

      id = mod (iday,30)
      im = int((iday - 1)/30.6) + 1
      enso_start = date2enso (id, im, iyear)

      call cprof(nx,ny,mask,i_prof,j_prof,ipt_prof)

!  new values, with ahmix0 included
      xk(1) = -65.218
      xk(2) = 1.8782
      xk(3) = -.00334
      ns   = inp_rarr('Spline_B'//CHAR(0), 3, xk,  xk)
      if (ns.ne.3) then
         print*,'need 3 values of Spline_B'
         stop
      endif
      b0 = xk(1)
      b1 = xk(2)
      b2 = xk(3)

      if (use_sst) then
         n_spline = inp_str('Spline_file'//CHAR(0),'DATA/spline-default.txt'//CHAR(0), fbspline)

         open (unit = 64, file = fbspline(1:n_spline))
         READ(64,*) nspk, nknots
         PRINT*, "nspk, nknots = ", nspk, nknots
         READ(64,*) knots(1:nknots)
         print*,'knots=',knots(1:nknots)
         READ(64,*) coeffs(1:nknots)
         print*,'coeffs=',coeffs(1:nknots)
         READ(64,*) HCUT, TCUT
         print*,'HCUT,TCUT=',HCUT,TCUT
         close (unit = 64)

         spknots = 1.d0*knots
         spcoeffs = 1.d0*coeffs
      endif

      return
      end

!----------------------------------------------------------------------
      subroutine make_iox(npt)
      use wrap
      implicit none
      include 'om_wrap.h'
!----------------------------------------------------------------------
      integer i, j, npt, nin

      nxyc = 0
      do i = 1, NX
      do j = 1, NY
         nin = (j-1)*NX + i
         if (mask(nin).ne.0) then
            nxyc = nxyc + 1
            iox(nxyc) = nin
         endif
      enddo
      enddo

      npt = nxyc 

      return
      end

!----------------------------------------------------------------------
      subroutine inout (ioin)
      use wrap
      implicit none
      include 'om_wrap.h'
!----------------------------------------------------------------------
!     unix routine for opening the input parameter file and the output
!     logfil from the command line arguments.


      integer ioin, lnblnk, num

      character*80 arg
      integer i, narg, iargc, in

      narg = iargc()

      i  = 0
      in = 0
      do while (i .lt. narg)
         i = i + 1
         call getarg (i, arg)
         if (arg .eq. '-h' .or. arg .eq. '-help') then
            goto 100
         elseif (arg .eq. '-i') then
            i = i + 1
            call getarg (i, finp)
            open (unit = ioin, file = finp)
            fout = finp(1:lnblnk(finp))//'.log'
            fcpu = finp(1:lnblnk(finp))//'.cpu'
            ftios = '.tios'//CHAR(0)
            in = 1
         elseif (arg .eq. '-o') then
            i = i + 1
            call getarg (i, fout)
         elseif (arg .eq. '-t') then
            i = i + 1
            call getarg (i, ftios)
         elseif (arg .eq. '-d') then
            i = i + 1
            call getarg (i, fbi)
            n_in = lnblnk(fbi)
            call dump_rstrt
            stop
         endif
      enddo

      if (in .eq. 0) goto 100

      open (unit = iout, file = fout)
      goto 200

  100 call getarg (0, arg)
      write (6, *) 
     *     'usage: '//arg(1:lnblnk(arg))//' [-i file][-t ftios][-d file]'
      write (6,*) 'where: -i file - for model control <file>'
      write (6,*) '       -t file - for tios control <file> (deflt:<.tios>)'
      write (6,*) '       -d file - make a dump of data/restart <file>'
      stop

  200 return
      end

      subroutine enso2res (renso, id, im, iy)
!----------------------------------------------------
#define LEAP_YEAR(y) (mod(y,4) .eq. 0)
      integer*2 norm(12)
      data norm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      integer int

      iy = int(renso/12.)

      res = renso - real(12*iy)
      im = int(res)

      res = abs(res - real(im))
      id = int(res*norm(im+1))

      return
      end


!----------------------------------------------------------------------
      function date2enso (id, im, iy)
      implicit none
!----------------------------------------------------------------------
#define LEAP_YEAR(y) (mod(y,4) .eq. 0)
      real date2enso
      integer id, im, iy
      integer*2 norm(12)
      data norm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      if (im .eq. 2 .and. LEAP_YEAR(iy) ) then
         date2enso = (iy - 1960.) * 12. + 1. + real((id-1))/29.
      else
         date2enso = (iy - 1960.) * 12. + real(im-1) + real((id-1))/norm(im)
      endif

      end

!----------------------------------------------------------------------
      subroutine enso2date (enso, id, im, iy)
      implicit none
!----------------------------------------------------------------------
#define LEAP_YEAR(y) (mod(y,4) .eq. 0)
      integer id,im,iy
      real enso, res

      integer*2 norm(12)
      data norm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      iy = 1960 + int(enso/12.)
      if (enso.lt.0) iy = iy - 1

      res = enso - (iy - 1960.)*12.
      im = int(res) + 1

      if (im.eq.13) then
         im = 1
         iy = iy + 1
      endif

      if (im .eq. 2 .and. LEAP_YEAR(iy) ) then
         id = 1 + int(29 * (res - int(res)))
      else
         id = 1 + int(norm(im) * (res - int(res)))
      endif

      end

      subroutine init_rstrt (npt,iseg)
!---------------------------------------------------------
      use wrap
      use wdata
      include 'om_wrap.h'
      include 'om_data.h'
      integer lnblnk
      character*100 str
      common /run/   nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast
      common /new_save/ iword, nruns, nscpu, nswll, inf(100), rnf(100)
      dimension   iseg(1)

      equivalence (inf(20),kd_sta),
     *            (inf(21),kd_xy),  (inf(23),kd_time), 
     *            (inf(24),kd_seg), (inf(25),kd_uv), (inf(26),kd_temp), 
     *            (inf(27),kd_salt),(inf(28),kd_tr), (inf(29),kd_phi), 
     *            (inf(30),kd_conv), (inf(31),kd_means),
     *            (inf(32),kd_psi), (inf(33),kd_ice),
     *            (inf(100),kd_end)

      call opda (1, 1, fbo(1:n_out))

      do i = 1, 100
         inf(i) = 0
         rnf(i) = 0.
      enddo

      iword   = MACHINE_WORD
      
      kd = iword * 200

      str = 'Restart file for '//finp(1:lnblnk(finp)) 
      call wrda (1, 1, kd, len(str), str)

      kd = iword * 200 + 100
      kd_sta = kd
      kd_xy  = kd
      call wrda (1, iword, kd, nx, xv) 
      call wrda (1, iword, kd, ny, yv) 

      call segm_from_iox (npt, nseg, iox, iseg, nx*ny)
      kd_seg = kd
      call wrda (1, iword, kd, 2*nseg, iseg) 
      
      kd_uv = kd

      nstart     = 1
      nlaststart = 1
      nlast      = nsteps
      nscpu = 0
      nswll = 0

      inf(1) = 120
      inf(2) = nx
      inf(3) = ny
      inf(4) = nmodes
      inf(5) = 1
      inf(6) = npt
      inf(7) = nseg

      inf(11) = nstart

      inf(15) = nlaststart
      inf(16) = nsteps
      inf(17) = nlast

      rnf(1) = delt
      rnf(2) = enso_start
      rnf(3) = enso_scale
      rnf(4) = rnf(2)

      return
      end

      subroutine read_rstrt (nxp, nyp, npt)
!--------------------------------------------------------------------
      use wrap
      include 'om_wrap.h'
      common /run/   nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast
      common /new_save/ iword, nruns, nscpu, nswll, inf(100), rnf(100)

      equivalence (inf(20),kd_sta),
     *            (inf(21),kd_xy),  (inf(23),kd_time), 
     *            (inf(24),kd_seg), (inf(25),kd_uv), (inf(26),kd_temp), 
     *            (inf(27),kd_salt),(inf(28),kd_tr), (inf(29),kd_phi), 
     *            (inf(30),kd_conv), (inf(31),kd_means),
     *            (inf(32),kd_psi), (inf(33),kd_ice),
     *            (inf(100),kd_end)

      iword = MACHINE_WORD
      
      call opda (1, 0, fbi(1:n_in))
      kd = 0
      call rdda (1, iword, kd, 100, inf) 
      kd = irdda (iword, 100, rnf) 

      ntype = idig(inf(1),1)
      if (ntype .ne. ITYPE_RST) then
         write (6, *) 'The file <' ,fbi(1:n_in),'> is not a RESTART-TYPE one !!'
         stop
      endif
      
      if (nxp.ne.inf(2) .or. nyp.ne.inf(3) .or.
     *    nmodes.ne.inf(4) .or. npt.ne.inf(6) ) then
         write(6,*)'The model dimensions are different from the restart data!'
         write(6,*)'nx=',nxp,inf(2),' ny=',nyp,inf(3),' nmodes=',nmodes,inf(4),' npt=',npt,inf(6)
         stop
      endif

      mmodess   = inf(6)*inf(4)

      nruns  = inf(12)

      if (irest.ne.4) then
         nstart = inf(11) + 1
         nscpu  = inf(13)
         nswll  = inf(14)
         if (irest .eq. 1) then
            nlaststart = inf(15)
         else
            nlaststart = nstart
         endif
         nlast = nlaststart + nsteps - 1
         enso_start = rnf(2)
         rnf(4) = rnf(2) + rnf(3)*nstart
      else
         nstart     = 1
         nlaststart = 1
         nlast      = nsteps
         nscpu = 0
         nswll = 0
         inf(11) = nstart

         inf(15) = nlaststart
         inf(16) = nsteps
         inf(17) = nlast

         rnf(2) = enso_start
         rnf(4) = rnf(2)
      endif

      if (irest .lt. 3) then
         ekf1 = rnf(10)
         epf1 = rnf(11)
         hcf1 = rnf(12)
         wcf1 = rnf(13)
         vlf1 = rnf(14)
      endif

      kd = kd_uv
      call rdda (1, iword, kd, mmodess, um) 
      kd = irdda (iword, mmodess, vm) 
      kd = irdda (iword, mmodess, vmy) 
      kd = irdda (iword, mmodess, hm) 
      kd = irdda (iword, nxy, sst)
      
      call clda(1)

      call opda (1, 1, fbo(1:n_out))
      inf(4)  = nmodes
      inf(12) = nruns + 1
      inf(15) = nlaststart
      inf(16) = nsteps
      inf(17) = nlast

      rnf(1) = delt

      return
      end

      subroutine keep_rstrt (nstep, nskip, tenso)
!--------------------------------------------------------------------
      use wrap
      include 'om_wrap.h'
      common /new_save/ iword, nruns, nscpu, nswll, inf(100), rnf(100)
      equivalence (inf(20),kd_sta),
     *            (inf(21),kd_xy),  (inf(23),kd_time), 
     *            (inf(24),kd_seg), (inf(25),kd_uv), (inf(26),kd_temp), 
     *            (inf(27),kd_salt),(inf(28),kd_tr), (inf(29),kd_phi), 
     *            (inf(30),kd_conv), (inf(31),kd_means),
     *            (inf(32),kd_psi), (inf(33),kd_ice),
     *            (inf(100),kd_end)

      if (nskip.eq.0) return
      if (mod(nstep, nskip) .ne. 0) goto 100

      print*,'saving restart file:',tenso

      npt   = inf(6)
      mmodess = npt*inf(4)

      kd = kd_uv
      call wrda (1, iword, kd, mmodess, um) 
      kd = iwrda (iword, mmodess, vm) 
      kd = iwrda (iword, mmodess, vmy) 
      kd = iwrda (iword, mmodess, hm) 
      kd = iwrda (iword, nxy, sst) 

      kd_end = kd

      inf(11) = nstep

      rnf(10) = ekf1 
      rnf(11) = epf1
      rnf(12) = hcf1
      rnf(13) = wcf1
      rnf(14) = vlf1

  100 continue
      inf(18) = nstep

      kd = 0
      call wrda (1, iword, kd, 100, inf) 
      kd = iwrda (iword, 100, rnf) 

      call flda (1)
      return
      end

      subroutine close_rstrt 
!---------------------------------
      call clda(1)
      return
      end

      subroutine dump_rstrt
!------------------------------------------------------------------
      use wrap
      include 'om_wrap.h'
      character*100 str
      character*4 mname(12)
      common /new_save/ iword, nruns, nscpu, nswll, inf(100), rnf(100)
      save mname

      data mname /'Jan.','Feb.','Mar.','Apr.','May ','Jun.',
     *            'Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'/

      iword = MACHINE_WORD
      
      call opda (1, 0, fbi(1:n_in))
      kd = 0
      call rdda (1, iword, kd, 100, inf) 
      kd = irdda (iword, 100, rnf)       
      kd = irdda (1, 100, str)       

      call clda(1)

      write (6, *) 'Dump of the data file <',fbi(1:n_in),'>:'
      write (6, *) 'LABEL: ', str
      write (6, *) 'TYPE', inf(1)
      write (6, *) 'NX =', inf(2)
      write (6, *) 'NY =', inf(3)
      write (6, *) 'NMODES =', inf(4)
      write (6, *) 'NT =', inf(5)
      write (6, *) 'NPACK =', inf(6)
      write (6, *) 'NSEGM =', inf(7)

      call enso2date (rnf(2), id, im, iy)
      write (6, 101) 'Model Starting Date:        ', mname(im), id, iy
      call enso2date (real(rnf(2)+rnf(3)*inf(17)), id, im, iy)
      write (6, 101) 'Scheduled End of the Run:   ', mname(im), id, iy
      call enso2date (real(rnf(2)+rnf(3)*inf(18)), id, im, iy)
      write (6, 101) 'Current End of the Run:   ' , mname(im), id, iy

      call enso2res (real(rnf(3)*inf(17)), id, im, iy)
      write (6, 201) 'Scheduled Length of the Run:', iy, im, id, inf(17)
      call enso2res (real(rnf(3)*inf(18)), id, im, iy)
      write (6, 201) 'Elapsed Model Time:         ', iy, im, id, inf(18)

      write (6, *) 'Number of restarts:', inf(12)

      call enso2date (rnf(4), id, im, iy)
      write (6, 101) 'Last Re-start:              ', mname(im), id, iy
      call enso2date (real(rnf(2) + rnf(3)*inf(15)), id, im, iy)
      write (6, 101) 'Last New Run:               ', mname(im), id, iy
      call enso2date (real(rnf(2) + rnf(3)*inf(11)), id, im, iy)
      write (6, 101) 'Last Save:                  ', mname(im), id, iy

      i = inf(13)
      if (i/3600 .eq. 1) then
         write (6, 301) 'Total CPU  time:', i/3600, mod(i,3600)/60, mod(i,60) 
      else
         write (6, 302) 'Total CPU  time:', i/3600, mod(i,3600)/60, mod(i,60) 
      endif
      i = inf(14)
      if (i/3600 .eq. 1) then
         write (6, 301) 'Total WALL time:', i/3600, mod(i,3600)/60, mod(i,60) 
      else
         write (6, 302) 'Total WALL time:', i/3600, mod(i,3600)/60, mod(i,60) 
      endif

  101 format (a30, a4, i2, ',' , i5)
  201 format (a30, i5, ' Years', i3, ' months', i3, ' days. (', i7, ' steps.)') 
  301 format (a20, i5, ' hour ', i3, ' min', i3, ' sec.')
  302 format (a20, i5, ' hours', i3, ' min', i3, ' sec.') 

      return
      end

      subroutine cprof (nx,ny,mask,i_prof,j_prof,ipt_prof)

      implicit none
      integer nx,ny,i_prof,j_prof,ipt_prof, i, j, mask(nx,ny), ipt

      ipt_prof = (j_prof-1)*nx + i_prof      

100   if (mask(i_prof,j_prof).eq.0) then
         print*,'nx,ny = ',nx,ny
         print*,'profile coordinates are not an ocean point', i_prof,j_prof
         stop
      endif

      return
      end

!----------------------------------------------------------------------
      subroutine read_mask (tag, nx, ny, mask)
!----------------------------------------------------------------------
      character*(*) tag
      character*1  ch, buff(119)
      character*64 number
      logical inp_def

      integer im,ip,jm,jp, mask(nx,ny), length

      equivalence (ch, buff(1)), (number, buff(2))

      if ( inp_def(tag) ) then
         ix = 1
         jy = ny
         ic0 = ichar('0')
         length = inp_wnxt(buff)
         do while ( length .gt. 0)
            if     (ch.eq.'0' .or. ch.eq.'1' .or. ch.eq.'2') then
               do i = 1, nx
                  mask(i,jy) = ichar(buff(i)) - ic0
               enddo
               jy = jy - 1
            elseif (ch.eq.'w' .or.ch.eq.'z' .or. ch.eq.'x') then
               read (number(1:length-1), *) kx 
               ix = ix + kx
               if (ix - 1 .gt. NX) goto 200
               if (ch.eq.'w') ich = 1 
               if (ch.eq.'x') ich = 2 
               if (ch.eq.'z') ich = 0 
               do i = ix-kx, ix-1
                  mask(i,jy) = ich
               enddo
               if (ix-1 .eq. NX) then
                  ix = 1
                  jy = jy - 1
               endif
            elseif (ch .eq. 'r') then
               read (number, *) ky
               jy = jy - ky + 2
               if (jy .lt. 1) goto 200
               do j = jy+ky-2, jy, -1
                  do i = 1, nx
                     mask(i,j) = mask(i,j+1)
                  enddo
               enddo
               jy = jy - 1
            else
               goto 200
            endif
            length = inp_wnxt(buff)
            if (jy .eq. 0) goto 100
         enddo

  100    continue

      else
         do j = 1, ny
            do i = 1, nx
               mask(i,j) = 1
            enddo
         enddo
      endif

      return

  200 write (6, *) '!!!read_mask: wrong mask data, i,j:', i,j
      stop
      end
