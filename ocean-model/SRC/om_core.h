      include 'om_para.h'

      integer  ival(MXY),jval(MXY)
     *        ,lrt_lc(MXY,MMODES), lbt_lc(MXY,MMODES), lot_lc(MXY,MMODES)
     *        ,lrt_rc(MXY,MMODES), lbt_rc(MXY,MMODES), lot_rc(MXY,MMODES)
     *        ,lrt_wb(MXY,MMODES), lbt_wb(MXY,MMODES), lot_wb(MXY,MMODES)
     *        ,lrt_eb(MXY,MMODES), lbt_eb(MXY,MMODES), lot_eb(MXY,MMODES)
     *        ,lptix(MXY), lptiy(MXY), lptixy(MXY), lpt_b(MXY)
     *        ,lpt_e(MXY), lpt_w(MXY), lpt_c(MXY), lpt_bc(MXY)
     *        ,lpt_n(MXY), lpt_s(MXY)
     *        ,lpt_ew(MXY),lpt_ns(MXY)
     *        ,it_l(NXMAX,NYMAX,MMODES), it_r(NXMAX,NYMAX,MMODES)
     *        ,l_xclose(MXY),l_yclose(MXY)
     *        ,lpt_eo(MXY), lpt_wo(MXY), lpt_no(MXY), lpt_so(MXY)
     *        ,lpt_neo(MXY), lpt_nwo(MXY), lpt_seo(MXY), lpt_swo(MXY)
      real  c_xclose(MXY),c_yclose(MXY)
      common /inc_indices/  ival, jval
     *                    , lrt_lc, lbt_lc, lot_lc
     *                    , lrt_rc, lbt_rc, lot_rc
     *                    , lrt_wb, lbt_wb, lrt_eb, lbt_eb, lpt_b, lpt_c, lpt_bc
     *                    , lptix, lptiy, lptixy, lpt_e, lpt_w, lpt_n, lpt_s
     *                    , lpt_ew, lpt_ns, it_l, it_r, l_xclose, l_yclose
     *                    , c_xclose, c_yclose
     *                    , lpt_eo, lpt_wo, lpt_no, lpt_so
     *                    , lpt_neo, lpt_nwo, lpt_seo, lpt_swo

      real at_l(NXMAX,NYMAX,MMODES), at_r(NXMAX,NYMAX,MMODES)
     *     ,dt_l(NXMAX,NYMAX,MMODES), dt_r(NXMAX,NYMAX,MMODES)
      common /inc_precompute/ at_l, at_r, dt_l, dt_r

      integer npt
     *                    , nrt_lc(MMODES), nbt_lc(MMODES), not_lc(MMODES)
     *                    , nrt_rc(MMODES), nbt_rc(MMODES), not_rc(MMODES)
     *                    , nrt_wb(MMODES), nbt_wb(MMODES), not_wb(MMODES)
     *                    , nrt_eb(MMODES), nbt_eb(MMODES), not_eb(MMODES)
     *                    , nptix, nptiy, nptixy, npt_e, npt_w, npt_n, npt_s
     *                    , npt_ew, npt_ns, n_xclose, n_yclose, npt_b
     *                    , npt_c, npt_bc, npt_eo, npt_wo, npt_no, npt_so
     *                    , npt_neo, npt_nwo, npt_seo, npt_swo
      common /inc_parameters/ npt
     *                    , nrt_lc, nbt_lc, not_lc
     *                    , nrt_rc, nbt_rc, not_rc
     *                    , nrt_wb, nbt_wb, nrt_eb, nbt_eb
     *                    , nptix, nptiy, nptixy, npt_e, npt_w, npt_n, npt_s
     *                    , npt_ew, npt_ns, n_xclose, n_yclose
     *                    , not_eb, not_wb, npt_b
     *                    , npt_c, npt_bc, npt_eo, npt_wo, npt_no, npt_so
     *                    , npt_neo, npt_nwo, npt_seo, npt_swo

      real cn(NYMAX,NXMAX),dn(NYMAX,NXMAX),rhs(NYMAX,NXMAX)
      common /inc_temps/ cn, dn, rhs

      real fcor(NYMAX+1,MMODES), gamm(NYMAX,NXMAX,MMODES)
     *    ,gamc(NYMAX,NXMAX,MMODES),gamp(NYMAX,NXMAX,MMODES)
     *    ,weight_x(NYMAX,NXMAX),weight_y(NYMAX,NXMAX)
      common /inc_coef/ fcor,gamm,gamc,gamp,weight_x,weight_y

      real tau(MMODES),taui(MMODES),thyi(MMODES),hysq(MMODES),thxi(MMODES)
      real eps(MMODES), hxsq(MMODES), vnu(MMODES), dt2(MMODES)
      common /inc_grid/ tau,taui,thyi,hysq,thxi,eps,hxsq,vnu, dt2

