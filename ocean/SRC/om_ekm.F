c-------------------------------------------------------------------------- 
      subroutine ekman(nx,ny,nmodes,F,G,rs,fcor,um,vm,uscl,dscl,
     *     hvmsf,hmixx,hstarp,dx,dy,u1,v1,u_ekman,v_ekman) 
      use wsst
      implicit none 
      include 'om_sst.h'
c-------------------------------------------------------------------------- 
c
c   INPUT: 
c     nx,ny 
c     nmodes      number of modes 
c     ii,jj       compressed to i,j indices  (from om_sst.h)
c     F, G        wind stresses 
c     rs          friction parameter (1/s)
c     fcor        dimensional Coriolis parameter 
c     um(k)
c     vm(k)       modal solution to SWEs 
c     uscl(k)     velocity scale 
c     dscl(k)	  scaling for vertical displacement (m)
c     hvmsf       surface values of eigenfunctions 
c     hmixx       depth of surface layer 
c     hstarp      1./(vertical sum of lighthill coefs)  (1/m)
c     dx,dy       grid spacings (m) 
c 
c  OUTPUT: 
c     u_ekman,v_ekman       ekman flow field
c     u1,v1                 total velocity in surface layer  

      integer nx,ny,nmodes
      real rs,hs,hstarp(nx,ny),dx,dy,F(ny,nx),G(ny,nx),fcor(ny),
     *     dscl(nmodes,nx,ny),um(ny,nx,nmodes),vm(ny,nx,nmodes),
     *     uscl(nmodes), hvmsf(nmodes,nx,ny)
      real u1(nx,ny),v1(nx,ny),u_ekman(nx,ny),v_ekman(nx,ny), hmixx(nx,ny)

      integer i,j,k,kk,mode
      real hstar,fco,denom,denomi,facm,dudx,dvdy,hvm
 
c   first compute the surface Ekman flow: 
 
      do k = 1, npt_u
         i = ii(k)  
         j = jj(k)  
         hstar = 1./hmixx(i,j) - hstarp(i,j)
c         print*,'hstar,hmixx(i,j),hstarp(i,j)',hstar,hmixx(i,j),hstarp(i,j)
c         stop
         if (hstar.lt.0.) then
            print*,'hstar is negative in ekman! '
            print*,i,j,hmixx(i,j),1./hstarp(i,j),hstar
            stop
         endif
         fco = fcor(j) 
         denom = rs**2 + fco**2 
         denomi = hstar/denom 
         
         u_ekman(i,j) = denomi* (rs* F(j,i) + fco* G(j,i)) 
         v_ekman(i,j) = denomi* (rs* G(j,i) - fco* F(j,i)) 
      enddo 
 
c   set no normal flow boundary conditions: 
      call bcset2(nx,ny,u_ekman,v_ekman)

c   compute total surface flow: 

      do k = 1, npt_u
         i = ii(k)  
         j = jj(k)  
         u1(i,j) = u_ekman(i,j)
         v1(i,j) = v_ekman(i,j)
         do mode = 1, nmodes 
            hvm  = hvmsf(mode,i,j)
            facm = dscl(mode,i,j)*uscl(mode)*hvm
            u1(i,j) = u1(i,j) + um(j,i,mode)*facm 
            v1(i,j) = v1(i,j) + vm(j,i,mode)*facm 
         enddo 
      enddo 
 
      return 
      end  

