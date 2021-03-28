      MODULE wrap

      integer, dimension(:), allocatable :: iox, ii, jj
      integer, dimension(:), allocatable :: nmixx

      integer, dimension(:), allocatable :: mask
      real, dimension(:), allocatable :: rmask

      real, dimension(:), allocatable :: tp
      real, dimension(:), allocatable :: F,G,Q,xv,yv,fc,xm,ym
      real, dimension(:), allocatable :: zin,sst,w1,tdeep,qnet,qdiag,sst_terms,qfp
      real, dimension(:), allocatable :: u_ekman,v_ekman,sstt,uq,u1,v1,zim
      real, dimension(:), allocatable :: cldt,slrt,ut,vt,wt,aht,att
      real, dimension(:), allocatable :: qflxt,dQdTt,dQdft
      real, dimension(:), allocatable :: rlh,sh,qlw,hmod,htherm,Fa,Qa,Ga

      real, dimension(:), allocatable :: anbv2, temp, salt

      real, dimension(:), allocatable :: p,u,v,t,um,vm,hm,vmy,un,vn,hn
      real, dimension(:), allocatable :: hvm, vvm, hcm

      real, dimension(:), allocatable :: hscl, dscl, uscl
      real, dimension(:), allocatable :: ahmix0,ahmixp,hvmsf,cph
      real, dimension(:), allocatable :: cl,sealvl,hmixx,hstarp,zcut

      END MODULE wrap
