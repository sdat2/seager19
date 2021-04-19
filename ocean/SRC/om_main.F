c-------------------------------------------------------------------------- 
      program om_inc   
      use wrap
      use wsst
      implicit none
      include 'om_wrap.h'
      include 'om_sst.h'
c-------------------------------------------------------------------------- 
c     multimode linear equatorial ocean model, using the INC scheme
c
c        authors:     N. Henderson(Naik), with contributions from: 
c                                   B. Blumenthal (multimode)
c                                   R. Seager     (AML)
c        references:        
c                    M. Israeli, Naik, N. and Cane, M.A., 1998. 
c             "The INC scheme for the shallow water equations"

c                    M.B. Blumenthal and Cane, M., 1989. "Accounting for 
c              parameter uncertainties in model verification: an illustration
c              with tropical sea surface temperature," 
c              J. Phys. Oceanogr.19, 815-830.

c                    R. Seager, Blumenthal, M.B. and Kushnir, Y., 1995.
c             "An advective atmospheric mixed layer model for ocean
c	       modeling purposes: Global simulation of surface heat fluxes",
c              J. Climate, 8, 1951-1964.

c                    R. Seager, Kushnir, Y. and Cane, M.A., 1995.
c             "On heat flux boundary conditions for ocean models"
c              J. Phys. Oceanogr., 25, 3219-3230.

c        description:              
c              
c-------------------------------------------------------------------------- 
      common /run/ nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast
      integer nstart,nlaststart,nskip,nsteps,nergy,nskipo,nlast,ioin
      data ioin,iout /1,2/

      common/winds/mtau,matau,tausc,atau,froude,itau_cos,drag_fac,mass
      integer mtau,matau,itau_cos,mass
      real tausc,atau,froude,drag_fac

      real tenso, dtt, dt_days,area
      integer nstep,iday_curr,iday_new
      integer id, im, iy, i, mode, kpt, kypt

      logical NEWRUN, non_stable

      integer nin,inp_modes,npt,ipt,in1,in2,in3

      real pin(MAXNZ)
      common /vertt/ pin

      real phasev(MMODES),scale_L
      real time_mul,run_time,save_time,rs,gamma

      real dtk(MMODES), hxk(MMODES), hyk(MMODES), dxm, dym, slat, fcor(NYMAX)
      real rayl(MMODES), rlap(MMODES)

      call inout (ioin)
      if (n_debug.ne.0) print*,'after inout'
      call model_input (tenso,inp_modes,phasev,time_mul,run_time,save_time,
     *                  gamma,rs,f1prime)
      if (n_debug.ne.0) print*,'after model_input'
      call model_memory 
      call make_iox(npt) 
      if (n_debug.ne.0) print*,'after make_iox'

      call logfile
      if (n_debug.ne.0) print*,'after logfile'


      if (use_profile) then
c                   read vertical density profile:
         if (nvert_prof.eq.0) then
            call vdens_init0(nx,ny,npt,MAXNZ,nin,pin,anbv2,temp,hsfc,hmixx,
     *                        z_cut,zcut) 
            if (n_debug.ne.0) print*,'after vdens_init0'
         else
            call vdens_init(nx,ny,npt,MAXNZ,nin,pin,temp,salt,anbv2,hmixx,nmixx,
     *                      tp,z_cut,zcut) 
            if (n_debug.ne.0) print*,'after vdens_init'
         endif
c                   compute vertical modes:
         do ipt = 1, nx*ny
            in1 = (ipt-1)*MAXNZ + 1         
            in2 = (ipt-1)*nmodes + 1
            in3 = (ipt-1)*nmodes*nz + 1
            call vmodessub(nin,pin,anbv2(in1),temp(in1),nmodes,zcut(ipt),
     *           fb_zis(1:n_zis),z_top,zin,nz,hstarp(ipt),
     *           cph(in2),cl(in2),sealvl(in2),
     *           hvm(in3),vvm(in3),hcm(in3),
     *           hvmsf(in2),ahmix0(ipt),ahmixp(in2))
         enddo
         call modepa(nx,ny,nmodes,npt,time_mul,delt,dxm,dym,xm,ym,X_MIN,X_MAX,
     *        Y_MIN,Y_MAX,cph,phasev,cl,sealvl,dtk,A_rayl,rayl,A_lapl,
     *        rlap,hxk,hyk,fcor,fc,hscl,uscl,dscl,ahmixp,scale_L)
      else
         call single_mode(nx,ny,time_mul,delt,dxm,dym,xm,ym,
     *        X_MIN,X_MAX,Y_MIN,Y_MAX,c_1mode,h_1mode,h1_1mode,
     *        hvm,hvmsf,dtk,A_rayl,rayl,hxk,hyk,fcor,fc,
     *        dscl,hscl,uscl,ahmix0,ahmixp,hstarp,area,scale_L)
         do ipt = 1, nx*ny
            hmixx(ipt) = hsfc
         enddo
 
      endif

      if (n_debug.ne.0) print*,'after mode initialization'

      dt_days = delt/3600./24.
      enso_scale = dt_days*12./365.     ! (enso_months/time_step)

      nsteps = run_time/delt
      nskip = save_time/delt


      nlast = nsteps
      if (irest .eq. 0) then
         NEWRUN = .true. 
         call init_rstrt(nxy,tp)
      else
         NEWRUN = irest .ge. 3
         if (irest.eq.1 .or. irest.eq.2) CALL TIOS_CNTRL (3, 1)
         call read_rstrt (nx, ny, nxy)
      endif

      if (n_debug.ne.0) print*,'after restart work'

      call init_data_out (ftios, fbt)
      if (n_debug.ne.0) print*,'after init_data_out'

      dtt = 12.*dt_days/365.
      call tios_cntrl (4, 10.*dtt)
      if (n_debug.ne.0) print*,'after tios_cntrl'


      iday_new  = int(nstart*dt_days)
      iday_curr = iday_new

      iflag = .false.

      do mode = 1, nmodes
         kpt = (mode-1)*NXY + 1
         kypt= (mode-1)*NY + 1
         call inc_init(nx,ny,npt,mode,mask,dtk(mode),hxk(mode),hyk(mode),
     *       fc(kypt),rayl(mode),rlap(mode),um(kpt),vm(kpt),hm(kpt),
     *       un(kpt),vn(kpt),hn(kpt),vmy(kpt),ischeme,ifper)
      enddo
      if (n_debug.ne.0) print*,'after inc_init'

      call forcing_init (nstart,nx,ny,xv,yv,area,F,G,Q,tp,hm,ym,scale_L)

      if (n_debug.ne.0) print*,'after forcing_init'
      call sst_index (nx,ny)
      if (n_debug.ne.0) print*,'after sst_index'

      if (use_sst) then
         call qflux_init (npt,nx,ny, xv,yv, iox,slat, sstt,
     *                      qflxt,dQdTt,dQdft,tp)
         if (n_debug.ne.0) print*,'after qflux_init' 
         call sst_init (nx,ny,sst,sstt,xv,yv,tp,irest)
         if (n_debug.ne.0) print*,'after sst_init' 
      endif

      tenso = enso_start + enso_scale * (nstart-1)

      do nstep = nstart, nlast
         
         tenso = enso_start + enso_scale * nstep
         
         if (mtau .eq. 1.or.mtau .eq. 6) then
            call forcing (nstep,nx,ny,npt,xv,yv,F,G,Q,tp,hm)
         endif

         if (n_debug.ne.0) print*,'after forcing' 

         do i = 1, nx*ny
            Fa(i) = F(i)
            Ga(i) = G(i)
            Qa(i) = Q(i)
         enddo

         do mode = 1, nmodes
            kpt = (mode-1)*NXY + 1
            if (ischeme.eq.1) then
               call inc(nstep,mode,nx,ny,hm(kpt),um(kpt),
     *              vm(kpt),vmy(kpt),Fa,Ga,Qa,ifper)
            elseif (ischeme.eq.2) then
               call leap_frog(nstep,mode,nx,ny,hm(kpt),um(kpt),vm(kpt),
     *                   hn(kpt),un(kpt),vn(kpt),Fa,Ga,Qa,ifper)
            endif
         enddo

         if (n_debug.ne.0) print*,'after inc time step' 

         if (use_sst) then
            call qflux(nstep, tenso, npt, nx, ny, iox, slat, qnet, 
     *       qfp,sst,sstt,qflxt,dQdTt,dQdft,f1prime,tp)
         endif
         
         if (use_ekm) then
            if (n_debug.ne.0) print*,'before ekman',nx,ny,nmodes,hmixx(1)

            call ekman(nx,ny,nmodes,F,G,rs,fcor,um,vm,uscl,dscl,
     *                 hvmsf,hmixx,hstarp,dxm,dym,u1,v1,u_ekman,v_ekman) 
            if (n_debug.ne.0) print*,'after ekman' 
            call upwelling(nx,ny,dxm,dym,hsfc,u1,v1,w1,ifper)
            if (n_debug.ne.0) print*,'after upwelling' 
         endif

         call tdeep_calc(nx,ny,nmodes,ahmix0,ahmixp,hm,tdeep,hmod,htherm)

         if (use_sst) then
            if (n_debug.ne.0) print*,'after tdeep_calc' 
            call sst_update(nx,ny,delt,dxm,dym,hsfc,gamma,u1,v1,w1,sst,sstt,
     *           qnet,qflxt,dQdTt,dQdft,tdeep,tp,tp(nxy+1),yv,idiag,qdiag,
     *           sst_terms,ifper,mask)

            if (n_debug.ne.0) print*,'after sst_update' 
         endif

         iday_new = int(nstep*dt_days)
         call keep_rstrt(nstep, nskip, tenso)

         if (iflag) then
            write (iout, *) 'Stable:ERROR, step', nstep
            goto 333
         endif
         
         if (NEWRUN .or. (iday_new .ne. iday_curr) ) then
            iday_curr = iday_new
            if (NEWRUN) NEWRUN = .false.
         endif
         
         call flush(iout)

         if (n_debug.ne.0) print*,'before data_out' 
         call data_out (tenso)     ! we also check stability
         if (n_debug.ne.0) print*,'after data_out' 
         
      enddo
      
      goto 444
c.............ABnormally finished run: 
  333 call tios_cntrl (1, 1)
      call data_out (tenso)

c.............normally finished run: 
  444 call close_rstrt
      write (iout, *) 'Finished at step =', nstep
      call enso2date (tenso, id, im, iy)
      write (iout, *) '<enso time> <', tenso, '>' 
      write (iout, *) '<Day:Month:Year> <',id,':',im,':',iy,'>' 

      print*, 'Finished at step =', nstep
      print*, '<enso time> <', tenso, '>' 
      print*, '<Day:Month:Year> <',id,':',im,':',iy,'>' 
      
      end

