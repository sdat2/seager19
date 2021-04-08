c---------------------------------------------
      subroutine qflux(nstep, tenso, npt, nx, ny, iox, slat, qnet,
     *                 qfp,sst, sstt, qflxt, dQdTt, dQdft, f1prime, tp)

      use wdata
      include 'om_data.h'
c---------------------------------------------
c.....AML(TCC - linearized heat fluxes from data)
c        sst     calculated sst
c        sstf    sst data
c        qflxf    qflx data
c        f1prime  change in f1 (greenhouse gas forcing) from 1958 to 2014

      integer iox(npt)
      real sst(1), qnet(1), sstt(1), ut(1),vt(1),wt(1), qflxt(1),
     *     dQdTt(1), dQdft(1), qfp(1),
     *     aht(1), att(1), cldt(1), slrt(1), sh(1), rlh(1), qlw(1), tp(1)
      real f1prime, tstart, tend, r

      call it_catch (nsst, tsst, nstep, it1, it2, tscl)

      NXY = nx*ny

      if (it2 .ne. isst) then
         isst = it2
         
         do i = 1, NXY
            sstf(i,1) = sstf(i,2)
         enddo
         
         nxy = nx*ny
         call read_zt (idf_sst, lsst, 1, it2, 'sst', tp, sstf(1,2)) 
      endif

      do k = 1, NXY
         sstf(k,3) = sstf(k,1) + tscl * (sstf(k,2) - sstf(k,1))
      enddo

      call it_catch (nqflx, tqflx, nstep, it1, it2, tscl)
      if (it2 .ne. iqflx) then
         iqflx = it2
         
         do i = 1, NXY
            qflxf(i,1) = qflxf(i,2)
         enddo
         
         nxy = nx*ny
         call read_zt (idf_qflx, lqflx, 1, it2, 'qflx', tp, qflxf(1,2)) 
      endif

      do k = 1, NXY
         qflxf(k,3) = qflxf(k,1) + tscl * (qflxf(k,2) - qflxf(k,1))
      enddo
         
      call it_catch (ndQdT, tdQdT, nstep, it1, it2, tscl)
      if (it2 .ne. idQdT) then
         idQdT = it2
         
         do i = 1, NXY
            dQdTf(i,1) = dQdTf(i,2)
         enddo
         
         nxy = nx*ny
         call read_zt (idf_dQdT, ldQdT, 1, it2, 'dQdT', tp, dQdTf(1,2)) 
      endif

      do k = 1, NXY
         dQdTf(k,3) = dQdTf(k,1) + tscl * (dQdTf(k,2) - dQdTf(k,1))
      enddo
         
      call it_catch (ndQdf, tdQdf, nstep, it1, it2, tscl)
      if (it2 .ne. idQdf) then
         idQdf = it2
         
         do i = 1, NXY
            dQdff(i,1) = dQdff(i,2)
         enddo
         
         nxy = nx*ny
         call read_zt (idf_dQdf, ldQdf, 1, it2, 'dQdf', tp, dQdff(1,2)) 
      endif

      do k = 1, NXY
         dQdff(k,3) = dQdff(k,1) + tscl * (dQdff(k,2) - dQdff(k,1))
      enddo
         
      call interp_u2t(nx,ny,sstf(1,3),sstt)
      call interp_u2t(nx,ny,qflxf(1,3),qflxt)
      call interp_u2t(nx,ny,dQdTf(1,3),dQdTt)
      call interp_u2t(nx,ny,dQdff(1,3),dQdft)

c   qnet = dQdT * (t1 - sst) + dQdf * fprime
c     where fprime = f1p * (t-tstart)/(tend-tstart)
c     and f1p ramps in time from 0 to f1prime (~= -0.006)

      tstart = date2enso(1,1,1958)
      tend   = date2enso(1,1,2015)
      r = max((tenso - tstart)/(tend - tstart),0.0)

      do k = 1, NXY
        qnet(k) = dQdTt(k)*(sst(k)-sstt(k))+ dQdft(k)*r*f1prime
        qfp(k) = dQdft(k)*r*f1prime
      enddo

      return
      end

c---------------------------------------------
      subroutine qflux_init (npt, nx, ny, xv, yv, iox, slat,
     *                      sstt,qflxt,dQdTt,dQdft,tp)
      use wdata
      include 'om_data.h'
c---------------------------------------------
      integer iox(npt), tp(1)
      real xv(nx),yv(ny),sstt(nx,ny)
      real qflxt(nx,ny), dQdTt(nx,ny), dQdft(nx,ny)
            
      nxy = nx*ny

c--------------------------------------------------
      allocate (sstf(nxy,3))
      call odb_open(idf_sst, fbsst(1:n_sst), 0)
      call odb_rddm(idf_sst, 'T', nsst)

      allocate (tsst(nsst))
      call odb_rdgr(idf_sst, 'T', nsst, tsst)

      call it_catch (nsst, tsst, nstart, it1, it2, tscl)
      isst = it2
      call data_on_model_grid(nx,ny,npt,idf_sst, lsst, 'sst', tp, 1)
      
      call read_zt (idf_sst, lsst, 1, it1, 'sst', tp, sstf(1,1)) 
      call read_zt (idf_sst, lsst, 1, it2, 'sst', tp, sstf(1,2)) 

      do k = 1, nxy
         sstf(k,3) = sstf(k,1) + tscl * (sstf(k,2) - sstf(k,1)) 
      enddo

c--------------------------------------------------
      allocate (qflxf(nxy,3))
      call odb_open(idf_qflx, fbqflx(1:n_qflx), 0)
      call odb_rddm(idf_qflx, 'T', nqflx)

      allocate (tqflx(nqflx))
      call odb_rdgr(idf_qflx, 'T', nqflx, tqflx)

      call it_catch (nqflx, tqflx, nstart, it1, it2, tscl)
      iqflx = it2
      call data_on_model_grid(nx,ny,npt,idf_qflx, lqflx, 'qflx', tp, 1)
      
      call read_zt (idf_qflx, lqflx, 1, it1, 'qflx', tp, qflxf(1,1)) 
      call read_zt (idf_qflx, lqflx, 1, it2, 'qflx', tp, qflxf(1,2)) 

      do k = 1, nxy
         qflxf(k,3) = qflxf(k,1) + tscl * (qflxf(k,2) - qflxf(k,1)) 
      enddo
      
c--------------------------------------------------
      allocate (dQdTf(nxy,3))
      call odb_open(idf_dQdT, fbdQdT(1:n_dQdT), 0)
      call odb_rddm(idf_dQdT, 'T', ndQdT)

      allocate (tdQdT(ndQdT))
      call odb_rdgr(idf_dQdT, 'T', ndQdT, tdQdT)

      call it_catch (ndQdT, tdQdT, nstart, it1, it2, tscl)
      idQdT = it2
      call data_on_model_grid(nx,ny,npt,idf_dQdT, ldQdT, 'dQdT', tp, 1)
      
      call read_zt (idf_dQdT, ldQdT, 1, it1, 'dQdT', tp, dQdTf(1,1)) 
      call read_zt (idf_dQdT, ldQdT, 1, it2, 'dQdT', tp, dQdTf(1,2)) 

      do k = 1, nxy
         dQdTf(k,3) = dQdTf(k,1) + tscl * (dQdTf(k,2) - dQdTf(k,1)) 
      enddo
      
c-------------------------------------------------
      allocate (dQdff(nxy,3))
      call odb_open(idf_dQdf, fbdQdf(1:n_dQdf), 0)
      call odb_rddm(idf_dQdf, 'T', ndQdf)

      allocate (tdQdf(ndQdf))
      call odb_rdgr(idf_dQdf, 'T', ndQdf, tdQdf)

      call it_catch (ndQdf, tdQdf, nstart, it1, it2, tscl)
      idQdf = it2
      call data_on_model_grid(nx,ny,npt,idf_dQdf, ldQdf, 'dQdf', tp, 1)
      
      call read_zt (idf_dQdf, ldQdf, 1, it1, 'dQdf', tp, dQdff(1,1)) 
      call read_zt (idf_dQdf, ldQdf, 1, it2, 'dQdf', tp, dQdff(1,2)) 

      do k = 1, nxy
         dQdff(k,3) = dQdff(k,1) + tscl * (dQdff(k,2) - dQdff(k,1)) 
      enddo
      
      nxm = nx
      nym = ny
      allocate (dxd(nxm))
      dxd = 0.
      allocate (dyd(nym))
      dyd = 0.

      do j = 1, nym
         dyd(j) = yv(2) - yv(1)
      enddo

      do i = 1, nxm
         dxd(i) = xv(2) - xv(1)
      enddo
      
      slat = yv(1) + dyd(1)/2.

      call interp_u2t(nx,ny,sstf(1,3),sstt)
      call interp_u2t(nx,ny,qflxf(1,3),qflxt)
      call interp_u2t(nx,ny,dQdTf(1,3),dQdTt)
      call interp_u2t(nx,ny,dQdff(1,3),dQdft)

      return
      end
      
c---------------------------------------------------------------------------
      subroutine interp_u2t(nx,ny,fu,ft)
      use wsst
      implicit none 
      include 'om_sst.h'
c---------------------------------------------------------------------------
      integer nx,ny
      real fu(nx,ny),ft(nx,ny) 

      integer i, j

      do i = 1, nx
      do j = 1, ny
         ft(i,j) = fu(i,j)
      enddo
      enddo

      return
      end
