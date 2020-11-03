      MODULE wdata

      real, dimension(:), allocatable :: taux0,tauy0,taux1,tauy1
      real, dimension(:), allocatable :: ttau, tsst, tprcp, tclm, zclm
      real, dimension(:), allocatable :: tqflx, tdQdT, tdQdf
      real, dimension(:), allocatable :: xd, yd, blcf
      real, dimension(:), allocatable :: dxd,dyd
      integer, dimension(:), allocatable :: im2d,ixd,itp
      real, dimension(:,:), allocatable :: sstf,slrf,cldf,uwnd,vwnd,wnsp
      real, dimension(:,:), allocatable :: qflxf, dQdTf, dQdff
      real, dimension(:,:), allocatable :: ahum,atem
      real, dimension(:), allocatable :: amhum,amth

      END MODULE wdata
