c-------------------------------------------------
      subroutine model_memory 
      use wrap
      implicit none
      include 'om_wrap.h'
c-------------------------------------------------
      integer nxyz,mxyz,mzm,mxyzp,i
 
      nxy  = nx*ny
      nxyz = nxy * nz
      mxyz = nxy * nmodes
      mxyzp = (nxy + nx) * nmodes
      mzm  = nz * nmodes

c  3D---------------------------------------------
      allocate (tp(nx*ny*MAXNZ))
      tp = 0.

      allocate (u(nxyz))
      u = 0.
      allocate (v(nxyz))
      v = 0.
      allocate (p(nxyz))
      p = 0.
      allocate (t(nxyz))
      t = 0.

      allocate (um(mxyz))
      um = 0.
      allocate (vm(mxyz))
      vm = 0.
      allocate (hm(mxyz))
      hm = 0.
      allocate (vmy(mxyz))
      vmy = 0.

      allocate (hn(mxyz))
      hn = 0.
      allocate (un(mxyz))
      un = 0.
      allocate (vn(mxyz))
      vn = 0.


c  2D---------------------------------------------

      allocate (tdeep(nxy))
      allocate (htherm(nxy))
      allocate (hmod(nxy))
      allocate (u1(nxy))
      u1 = 0.
      allocate (v1(nxy))
      v1 = 0.
      allocate (u_ekman(nxy))
      allocate (v_ekman(nxy))
      allocate (sst(nxy))
      allocate (sst_terms(7*nxy))
      allocate (w1(nxy))
      allocate (qnet(nxy))
      allocate (qfp(nxy))
      allocate (qdiag(nxy))
      allocate (sstt(nxy))
      allocate (qflxt(nxy))
      allocate (dQdTt(nxy))
      allocate (dQdft(nxy))

      allocate (cldt(nxy))
      allocate (slrt(nxy))
      allocate (ut(nxy))
      allocate (vt(nxy))
      allocate (wt(nxy))
      allocate (aht(nxy))
      allocate (att(nxy))
      allocate (rlh(nxy))
      allocate (sh(nxy))
      allocate (qlw(nxy))

      allocate (F(nxy))
      allocate (G(nxy))
      allocate (Q(nxy))
      allocate (Fa(nxy))
      allocate (Ga(nxy))
      allocate (Qa(nxy))

      allocate (xm(nxy))
      allocate (ym(nxy))

      allocate (fc(ny*nmodes))

      allocate (iox(nxy))

c  1D---------------------------------------------

      allocate (hvm(nz*nmodes*nxy))
      allocate (vvm(nz*nmodes*nxy))
      allocate (hcm(nz*nmodes*nxy))

      allocate (anbv2(MAXNZ*nxy))
      allocate (temp(MAXNZ*nxy))
      allocate (salt(MAXNZ*nxy))

      allocate (hvmsf(nxy*nmodes))
      allocate (ahmixp(nxy*nmodes))
      allocate (ahmix0(nxy))
      allocate (dscl(nxy*nmodes))
      allocate (uscl(nmodes))
      allocate (hscl(nmodes))

      allocate (cph(nxy*nmodes))
      allocate (cl(nxy*nmodes))
      allocate (sealvl(nxy*nmodes))

      allocate (hmixx(nxy))
      allocate (hstarp(nxy))
      allocate (nmixx(nxy))

      allocate (zcut(nxy))

      do i = 1, nxy
         nmixx(i) = 1
      enddo

      return
      end

c--------------------------------------
      subroutine datagrid_memory (nx,ny,npt,tmp)
      use wdata
      include 'om_data.h'
c--------------------------------------
      integer npt
      dimension tmp(1)
 
      allocate (xd(mxp))
      do i = 1, mxp
         xd(i) = tmp(i)
      enddo

      allocate (yd(myp))
      do i = 1, myp
         yd(i) = tmp(mxp+i)
      enddo

      npt2 = npt + 1

      allocate (ixd(mxp*myp))
      allocate (im2d(npt2))
      allocate (blcf(4*npt2))

c      print*,'datagrid_memory:',mxp,xd(1),xd(mxp)
c      print*,'datagrid_memory:',myp,yd(1),yd(myp)

      idatgr = 1
      
      return
      end
