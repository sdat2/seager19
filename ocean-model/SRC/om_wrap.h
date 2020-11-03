      include 'om_para.h'

      real enso_start, enso_scale, delt, delx, dely, ymin, hsfc
      common /MOD_CONST/ enso_start, enso_scale, delt, delx, dely, ymin, hsfc

      integer iout, n_in, n_out, n_tios, n_dir, n_zis,nz_mixed, n_debug
      character*80 ftios, fout, fcpu, finp, fbt, fbo, fbi, fbdir, fb_zis
      common /MOD_FILES/ ftios, fout, fcpu, finp, fbdir, iout
     *        , n_in, n_out, n_tios, fbt, fbo, fbi, n_dir, n_debug
     *        , fb_zis,n_zis,nz_mixed

      logical ifper
      integer NX, NY, NXY, nbx, nz, nbc, nmodes, nxyc, nxyct
      common /MOD_GRID/ NX, NY, NXY, nbx, nz, nbc, nmodes, nxyc, nxyct
      common /MOD_GRID3/ ifper

      real X_MAX, Y_MAX, X_MIN, Y_MIN
      common /MOD_GRID2/ X_MAX, Y_MAX, X_MIN, Y_MIN

      integer irest, ischeme, idiag
      common /new_misc/  irest, ischeme, idiag

      real c_1mode, H_1mode, H1_1mode, A_rayl, A_lapl 
      common /SEA_param/ c_1mode, H_1mode, H1_1mode, A_rayl,A_lapl

      logical use_profile
      integer nvert_prof
      real xlon_prof, ylat_prof, hmix_param
      integer ipt_prof
      common /vert_profile/ use_profile,xlon_prof, ylat_prof, nvert_prof,
     *                      ipt_prof, hmix_param

      real z_top,z_cut,z_bot,dzm

      common /vert_mode/z_top,z_cut,z_bot,dzm

      logical iflag
      common /stable/iflag

      real trans_coef_sst,rlx_time_sst 
      integer msst, mekm
      logical use_sst, use_ekm
      common /hflx_forc/msst,use_sst,use_ekm
     *                 ,trans_coef_sst,rlx_time_sst
