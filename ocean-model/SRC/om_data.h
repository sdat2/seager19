      character*80  fbden, fbwnd,fbtem,fbsst,fbprp,fbq,fbcld,fbslr,fbpot,fbsal
     *             ,fbqflx, fbdQdT, fbdQdf, fbmask, fbspline
      integer  n_den,n_wnd,n_tem,n_sst,n_slr,n_prp,n_q,n_cld,n_pot,n_sal
     *        , n_qflx, n_dQdT, n_dQdf, n_mask, n_spline
      common /new_files/ n_den, n_wnd,n_tem,n_sst, n_slr, n_prp, n_q, n_cld, 
     *                   n_pot, n_sal,
     *                   fbden, fbwnd,fbtem,fbsst,fbprp,fbq,fbcld,fbslr,
     *                   fbpot, fbsal
     *                  ,fbqflx, n_qflx
     *                  ,fbdQdT, n_dQdT
     *                  ,fbdQdf, n_dQdf, n_mask, fbmask, n_spline, fbspline
 
      integer  idf_dp, idf_cld, idf_slr, idf_sst,
     *     idf_tx, idf_ty,  ltau,itau,ntau,
     *     lsst,isst,nsst, idf_prp, lprp,iprp,nprp, idf_q, lq, iq, nq, 
     *     idf_t,idf_s,idf_hcl, lclm,iclm, ntclm, nzclm, lahum
     *    ,idf_qflx,lqflx,iqflx,nqflx
     *    ,idf_dQdT,ldQdT,idQdT,ndQdT
     *    ,idf_dQdf,ldQdf,idQdf,ndQdf, idf_mask

      real  cld_tscl, slr_tscl, tau_tscl, sst_tscl, clm_tscl

      common /new_forc/  idf_dp, idf_cld, idf_slr, idf_sst, cld_tscl, slr_tscl,
     *     idf_tx, idf_ty,  ltau,itau,ntau,tau_tscl, 
     *     lsst,isst,nsst, sst_tscl, lahum,
     *     idf_prp, lprp,iprp,nprp, idf_q, lq, iq, nq,
     *     idf_t,idf_s,idf_hcl, lclm,iclm, ntclm, nzclm,clm_tscl
     *    ,idf_qflx,lqflx,iqflx,nqflx
     *    ,idf_dQdT,ldQdT,idQdT,ndQdT
     *    ,idf_dQdf,ldQdf,idQdf,ndQdf, idf_mask
     

      integer  n_wsp, n_uwd, n_vwd, n_ah, n_at, n_prec,
     *         idf_wsp, idf_uwd, idf_vwd, idf_ah, idf_at, idf_prec
      character*80       fwsp, fuwd, fvwd, fah, fat, fprec
      common /pbl_files/ n_wsp, n_uwd, n_vwd, n_ah, n_at, n_prec,
     *                   idf_wsp, idf_uwd, idf_vwd, idf_ah, idf_at, idf_prec,
     *                   fwsp, fuwd, fvwd, fah, fat, fprec


      integer idatgr, mpack, mseg, mxp, myp, msx, msy
      common /forcgr/ idatgr, mpack,mseg, mxp,myp, msx,msy
 

      integer i_east, j_east
      real enso_k, enso_mu, enso_x1_fac, enso_x2_fac, enso_n
      real enso_x1, enso_x2
      common /enso/ enso_k, enso_mu, enso_x1_fac, enso_x2_fac, enso_n
     * , i_east, j_east, enso_x1, enso_x2
