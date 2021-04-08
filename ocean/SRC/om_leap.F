c----------------------------------------------------------------------
      subroutine leap_frog (nstep,k,nxr,nyr,hm,um,vm,hn,un,vn,F,G,Q,ifper)
      implicit none
      include 'om_core.h'
c----------------------------------------------------------------------
c
c     time steps the system:
c
c        (up-um)/2dt - fcor*vn  + d(hn)/dx = F - r*un + nu*(un_xx + un_yy)
c
c        (vp-vm)/2dt + fcor*un  + d(hn)/dy = G - r*vn + nu*(vn_xx + vn_yy)
c
c        (hp-hm)/2dt + d(un)/dx + d(vn)/dy = Q - r*hn + nu*(hn_xx + hn_yy)
c
c      on a C-grid
c
c
c     input:
c        nxr      : number of grid points in zonal direction
c        nyr      : number of grid points in meridional direction
c        hm,um,vm : k-th scaled mode, at t = (n-1)*dt
c        hn,un,vn : k-th scaled mode, at t = n*dt
c        F,G      : momentum forcing
c        Q        : mass forcing

      integer nstep, k, nxr, nyr
      logical ifper
      real F(nyr,nxr),G(nyr,nxr),Q(nyr,nxr)
      real hm(nyr,nxr),um(nyr,nxr),vm(nyr,nxr)
      real hn(nyr,nxr),un(nyr,nxr),vn(nyr,nxr)

      integer i,j,ipt,ik
      real dtt2, hxi, hyi, fc, coefx, coefy, hx2, hy2, fm
      real uc(NYMAX,NXMAX),vc(NYMAX,NXMAX),u,v

      dtt2 = dt2(k)
      hxi = 2.*thxi(k)
      hyi = 2.*thyi(k)
      coefx = vnu(k)
      coefy = vnu(k)/10.
      hx2 = 1./hxsq(k)
      hy2 = 1./hysq(k)

      do ik = 1, nptixy
         ipt = lptixy(ik)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)     
         fm = 0.5*(fcor(j,k)+fcor(j-1,k))
         v = 0.25*(vn(j,i)+vn(j+1,i)+vn(j,i-1)+vn(j+1,i-1))
         u = 0.25*(un(j,i)+un(j,i+1)+un(j-1,i)+un(j-1,i+1))
         F(j,i) = F(j,i) + fc* v
         G(j,i) = G(j,i) - fm* u

         F(j,i) = F(j,i) - hxi*(hn(j,i) - hn(j,i-1))
         Q(j,i) = Q(j,i) - hxi*(un(j,i+1) - un(j,i))

         G(j,i) = G(j,i) - hyi*(hn(j,i) - hn(j-1,i))
         Q(j,i) = Q(j,i) - hyi*(vn(j+1,i) - vn(j,i))


         F(j,i) = F(j,i) + coefx*(um(j,i+1)+ um(j,i-1)- 2.*um(j,i))*hx2
         G(j,i) = G(j,i) + coefx*(vm(j,i+1)+ vm(j,i-1)- 2.*vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)+ hm(j,i-1)- 2.*hm(j,i))*hx2

         F(j,i) = F(j,i) + coefy*(um(j+1,i)+ um(j-1,i)- 2.*um(j,i))*hy2
         G(j,i) = G(j,i) + coefy*(vm(j+1,i)+ vm(j-1,i)- 2.*vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j+1,i)+ hm(j-1,i)- 2.*hm(j,i))*hy2
      enddo

      do ik = 1, npt_no
         ipt = lpt_no(ik)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)     
         fm = 0.5*(fcor(j,k)+fcor(j-1,k))
         v = 0.25*(vn(j,i)+vn(j,i-1))
         u = 0.25*(un(j,i)+un(j,i+1)+un(j-1,i)+un(j-1,i+1))
         F(j,i) = F(j,i) + fc* v
         G(j,i) = G(j,i) - fm* u

         F(j,i) = F(j,i) - hxi*(hn(j,i) - hn(j,i-1))
         Q(j,i) = Q(j,i) - hxi*(un(j,i+1) - un(j,i))

         G(j,i) = G(j,i) - hyi*(hn(j,i) - hn(j-1,i))
         Q(j,i) = Q(j,i) - hyi*(0.      - vn(j,i))

         F(j,i) = F(j,i) + coefx*(um(j,i+1)+ um(j,i-1)- 2.*um(j,i))*hx2
         G(j,i) = G(j,i) + coefx*(vm(j,i+1)+ vm(j,i-1)- 2.*vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)+ hm(j,i-1)- 2.*hm(j,i))*hx2
         F(j,i) = F(j,i) + coefy*(um(j-1,i)- um(j,i))*hy2
         G(j,i) = G(j,i) + coefy*(vm(j-1,i)- 2.*vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j-1,i)- hm(j,i))*hy2
      enddo
      do ik = 1, npt_so
         ipt = lpt_so(ik)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)     
         v = 0.25*(vn(j+1,i)+vn(j+1,i-1))
         F(j,i) = F(j,i) + fc* v

         F(j,i) = F(j,i) - hxi*(hn(j,i) - hn(j,i-1))
         Q(j,i) = Q(j,i) - hxi*(un(j,i+1) - un(j,i))

         Q(j,i) = Q(j,i) - hyi*(vn(j+1,i) - 0.)

         F(j,i) = F(j,i) + coefx*(um(j,i+1)+ um(j,i-1)- 2.*um(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)+ hm(j,i-1)- 2.*hm(j,i))*hx2
         F(j,i) = F(j,i) + coefy*(um(j+1,i)- um(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j+1,i)- hm(j,i))*hy2
      enddo
      do ik = 1, npt_eo
         ipt = lpt_eo(ik)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)     
         fm = 0.5*(fcor(j,k)+fcor(j-1,k))
         v = 0.25*(vn(j,i)+vn(j+1,i)+vn(j,i-1)+vn(j+1,i-1))
         u = 0.25*(un(j,i)+un(j-1,i))
         F(j,i) = F(j,i) + fc* v
         G(j,i) = G(j,i) - fm* u

         F(j,i) = F(j,i) - hxi*(hn(j,i) - hn(j,i-1))
         Q(j,i) = Q(j,i) - hxi*(0.      - un(j,i))

         G(j,i) = G(j,i) - hyi*(hn(j,i) - hn(j-1,i))
         Q(j,i) = Q(j,i) - hyi*(vn(j+1,i) - vn(j,i))

         F(j,i) = F(j,i) + coefx*(um(j,i-1)- 2.*um(j,i))*hx2
         G(j,i) = G(j,i) + coefx*(vm(j,i-1)- vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i-1)- hm(j,i))*hx2
         F(j,i) = F(j,i) + coefy*(um(j+1,i)+ um(j-1,i)- 2.*um(j,i))*hy2
         G(j,i) = G(j,i) + coefy*(vm(j+1,i)+ vm(j-1,i)- 2.*vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j+1,i)+ hm(j-1,i)- 2.*hm(j,i))*hy2
      enddo
      do ik = 1, npt_wo
         ipt = lpt_wo(ik)
         i = ival(ipt)
         j = jval(ipt)
         fm = 0.5*(fcor(j,k)+fcor(j-1,k))
         u = 0.25*(un(j,i+1)+un(j-1,i+1))
         G(j,i) = G(j,i) - fm* u

         Q(j,i) = Q(j,i) - hxi*(un(j,i+1) - 0.)

         G(j,i) = G(j,i) - hyi*(hn(j,i) - hn(j-1,i))
         Q(j,i) = Q(j,i) - hyi*(vn(j+1,i) - vn(j,i))

         G(j,i) = G(j,i) + coefx*(vm(j,i+1)- vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)- hm(j,i))*hx2
         G(j,i) = G(j,i) + coefy*(vm(j+1,i)+ vm(j-1,i)- 2.*vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j+1,i)+ hm(j-1,i)- 2.*hm(j,i))*hy2
      enddo
      do ik = 1, npt_neo
         ipt = lpt_neo(ik)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)     
         fm = 0.5*(fcor(j,k)+fcor(j-1,k))
         v = 0.25*(vn(j,i)+vn(j,i-1))
         u = 0.25*(un(j,i)+un(j-1,i))
         F(j,i) = F(j,i) + fc* v
         G(j,i) = G(j,i) - fm* u

         F(j,i) = F(j,i) - hxi*(hn(j,i) - hn(j,i-1))
         Q(j,i) = Q(j,i) - hxi*(0.      - un(j,i))

         G(j,i) = G(j,i) - hyi*(hn(j,i) - hn(j-1,i))
         Q(j,i) = Q(j,i) - hyi*(0.      - vn(j,i))

         F(j,i) = F(j,i) + coefx*(um(j,i-1)- 2.*um(j,i))*hx2
         G(j,i) = G(j,i) + coefx*(vm(j,i-1)- vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i-1)- hm(j,i))*hx2
         F(j,i) = F(j,i) + coefy*(um(j-1,i)- um(j,i))*hy2
         G(j,i) = G(j,i) + coefy*(vm(j-1,i)- 2.*vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j-1,i)- hm(j,i))*hy2
      enddo
      do ik = 1, npt_nwo
         ipt = lpt_nwo(ik)
         i = ival(ipt)
         j = jval(ipt)
         fm = 0.5*(fcor(j,k)+fcor(j-1,k))
         u = 0.25*(un(j,i+1)+un(j-1,i+1))
         G(j,i) = G(j,i) - fm* u

         Q(j,i) = Q(j,i) - hxi*(un(j,i+1) - 0.)

         G(j,i) = G(j,i) - hyi*(hn(j,i) - hn(j-1,i))
         Q(j,i) = Q(j,i) - hyi*(0.      - vn(j,i))

         G(j,i) = G(j,i) + coefx*(vm(j,i+1)- vm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)- hm(j,i))*hx2
         G(j,i) = G(j,i) + coefy*(vm(j-1,i)- 2.*vm(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j-1,i)- hm(j,i))*hy2
      enddo
      do ik = 1, npt_seo
         ipt = lpt_seo(ik)
         i = ival(ipt)
         j = jval(ipt)
         fc = fcor(j,k)     
         v = 0.25*(vn(j+1,i)+vn(j+1,i-1))
         F(j,i) = F(j,i) + fc* v

         F(j,i) = F(j,i) - hxi*(hn(j,i) - hn(j,i-1))
         Q(j,i) = Q(j,i) - hxi*(0.      - un(j,i))

         Q(j,i) = Q(j,i) - hyi*(vn(j+1,i) - 0.)

         F(j,i) = F(j,i) + coefx*(um(j,i-1)- 2.*um(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j,i-1)- hm(j,i))*hx2
         F(j,i) = F(j,i) + coefy*(um(j+1,i)- um(j,i))*hy2
         Q(j,i) = Q(j,i) + coefy*(hm(j+1,i)- hm(j,i))*hy2
      enddo
      do ik = 1, npt_swo
         ipt = lpt_swo(ik)
         i = ival(ipt)
         j = jval(ipt)

         Q(j,i) = Q(j,i) - hxi*(un(j,i+1) - 0.)
         Q(j,i) = Q(j,i) - hyi*(vn(j+1,i) - 0.)

         Q(j,i) = Q(j,i) + coefx*(hm(j,i+1)- hm(j,i))*hx2
         Q(j,i) = Q(j,i) + coefx*(hm(j+1,i)- hm(j,i))*hy2
      enddo

c  compute new time step (hp,up,vp) and store in (hn,un,vn)
      do ipt = 1, npt
         i = ival(ipt)
         j = jval(ipt)
         F(j,i) = um(j,i) + dtt2*F(j,i)
         G(j,i) = vm(j,i) + dtt2*G(j,i)
         Q(j,i) = hm(j,i) + dtt2*Q(j,i)
         um(j,i) = un(j,i)
         vm(j,i) = vn(j,i)
         hm(j,i) = hn(j,i)
         un(j,i) = F(j,i)
         vn(j,i) = G(j,i)
         hn(j,i) = Q(j,i)
      enddo

      call Rayleigh_friction(nxr,nyr,un,vn,hn,2.*eps(k))

      do ik = 1, npt_w
         ipt = lpt_w(ik)
         i = ival(ipt)
         j = jval(ipt)
         un(j,i) = 0.
      enddo
      do ik = 1, npt_s
         ipt = lpt_s(ik)
         i = ival(ipt)
         j = jval(ipt)
         vn(j,i) = 0.
      enddo

      if (ifper) then
         do j = 1, nyr
            hn(j,1)   = hn(j,nxr-1)
            hn(j,nxr) = hn(j,2)
            un(j,1)   = un(j,nxr-1)
            un(j,nxr) = un(j,2)
            vn(j,1)   = vn(j,nxr-1)
            vn(j,nxr) = vn(j,2)
            hm(j,1)   = hm(j,nxr-1)
            hm(j,nxr) = hm(j,2)
            um(j,1)   = um(j,nxr-1)
            um(j,nxr) = um(j,2)
            vm(j,1)   = vm(j,nxr-1)
            vm(j,nxr) = vm(j,2)
         enddo
      endif

      return
      end
      
