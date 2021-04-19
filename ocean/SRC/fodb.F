c---------------------------------------------------------------------------
c ODB is a FORTRAN-based simplified wrapping for NetCDF and HDF libraries, 
c
c It assumes that for HDF compatibility you should obtain a copy of "libcdfdf"
c library and compile with -DCDF_DF "cpp" option in order to create HDF.
c 
c Documentation available on:
c  http://rainbow.ldgo.columbia.edu/documentation/libraries/odb/readme.html
c Please report all bugs and problems to: senya@rainbow.ldeo.columbia.edu 
c
c Senya Basin, 1993-94. 

      subroutine odb_open (idf, file, key) 
c3----------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      integer idf0,idf
      character*(*) file
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr
      common /odb_loc2/ idf0,idfmode(100)

      idf0 = 65536 

      if (key .eq. 0) then
         idf = ncopn (file, NCNOWRIT, iok)
      elseif (key .eq. 1) then
         call odb_creat (file, len(file), 0)
         idf = ncopn (file, NCWRITE, iok)
         OMODE = NCSFIL (idf, NCNOFILL, RCODE)
      elseif (key .eq. 2) then
         call odb_creat (file, len(file), 1)
         idf = ncopn (file, NCWRITE, iok)
         OMODE = NCSFIL (idf, NCNOFILL, RCODE)
      endif


      call ncpopt (0)
      call ncredf (idf, iok)
      call ncpopt (NCVERBOS+NCFATAL)
      idfmode(idf/idf0+1) = 0

      IDTYPE = NCFLOAT

      return
      end

      subroutine odb_close (idf) 
c4------------------------------
      integer idf
      call ncclos(idf, iok)

      return
      end

      subroutine odb_dfend (idf)
c5-----------------------------------
      common /odb_loc2/ idf0,idfmode(100)

      if (idfmode(idf/idf0+1) .eq. 0) then
         call ncendf (idf, iok)
         idfmode(idf/idf0+1) = 1
      endif

      return
      end

      subroutine odb_dfsta (idf)
c6-----------------------------------
      common /odb_loc2/ idf0,idfmode(100)


      if (idfmode(idf/idf0+1) .eq. 1) then
         call ncsnc(idf, iok)
         call ncredf (idf, iok)
         idfmode(idf/idf0+1) = 0
      endif

      return
      end

      subroutine odb_dfdm (idf, dname, nn)
c7-------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) dname

      call odb_dfsta (idf)

      if (nn .gt. 0) then
         idd = ncddef (idf, dname, nn, iok)
      else
         idd = ncddef (idf, dname, NCUNLIM, iok)
      endif

      return
      end

      subroutine odb_dfgr (idf, gname, nn)
c8-------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) gname
      dimension idd(1)
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfsta (idf)

      call ncpopt (0)
      idv_curr = ncvid (idf, gname, iok) 
      call ncpopt (NCVERBOS+NCFATAL)

      if (iok .ne. 0 .or. idv_curr .eq. -1) then
         if (nn .gt. 0) then
            idd(1)   = ncddef (idf, gname, nn, iok)
         else
            idd(1)   = ncddef (idf, gname, NCUNLIM, iok)
         endif
         idv_curr = ncvdef (idf, gname, IDTYPE, 1, idd, iok) 
         idf_curr = idf
      endif

      return
      end

      subroutine odb_dftm (idf, gn1)
c9-----------------------------------
      character*(*) gn1

      call odb_dfgr (idf, gn1, 0)

      return
      end

      subroutine odb_dfvar1 (idf, gn1, vname)
c10---------------------------------------------
      character*(*) gn1, vname
      dimension idd(1)
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfsta (idf)

      idd(1) = ncdid (idf, gn1, iok)
      idv_curr = ncvdef (idf, vname, IDTYPE, 1, idd, iok) 
      idf_curr = idf

      return
      end

      subroutine odb_dfvar2 (idf, gn1, gn2, vname)
c11-------------------------------------------------
      character*(*) gn1, gn2, vname
      dimension idd(2)
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfsta (idf)

      idd(1) = ncdid (idf, gn1, iok)
      idd(2) = ncdid (idf, gn2, iok)

      idv_curr = ncvdef (idf, vname, IDTYPE, 2, idd, iok) 
      idf_curr = idf

      return
      end

      subroutine odb_dfvar3 (idf, gn1, gn2, gn3, vname)
c12-------------------------------------------------------
      character*(*) gn1, gn2, gn3, vname
      dimension idd(3)
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfsta (idf)

      idd(1) = ncdid (idf, gn1, iok)
      idd(2) = ncdid (idf, gn2, iok)
      idd(3) = ncdid (idf, gn3, iok)

      idv_curr = ncvdef (idf, vname, IDTYPE, 3, idd, iok) 
      idf_curr = idf

      return
      end

      subroutine odb_dfvar4 (idf, gn1, gn2, gn3, gn4, vname)
c13-----------------------------------------------------------
      character*(*) gn1, gn2, gn3, gn4, vname
      dimension idd(4)
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfsta (idf)

      idd(1) = ncdid (idf, gn1, iok)
      idd(2) = ncdid (idf, gn2, iok)
      idd(3) = ncdid (idf, gn3, iok)
      idd(4) = ncdid (idf, gn4, iok)

      idv_curr = ncvdef (idf, vname, IDTYPE, 4, idd, iok) 
      idf_curr = idf

      return
      end

      subroutine odb_dfvar (idf, ngr, vgrn, vname)
c14-------------------------------------------------
      character*(*) vgrn(1), vname
      dimension idd(100)
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfsta (idf)

      do i = 1, ngr
         idd(i) = ncdid (idf, vgrn(i), iok)
      enddo

      idv_curr = ncvdef (idf, vname, IDTYPE, ngr, idd, iok) 
      idf_curr = idf

      return
      end

      subroutine odb_wrgr (idf, gn1, xx)
c15-----------------------------------------
      character*(*) gn1
      dimension xx(1), ista(1), icnt(1)
      character*80 buff
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfend (idf)

      idd = ncdid (idf, gn1, iok)
      call ncdinq (idf, idd, buff, icnt, iok)
      
      idv_curr = ncvid (idf, gn1, iok) 
      ista(1) = 1
      call ncvpt (idf, idv_curr, ista, icnt, xx, iok) 
      call ncsnc (idf, iok)

      idf_curr = idf

      return
      end

      subroutine odb_wrtm (idf, tname, indx, tval)
c16-------------------------------------------------
      character*(*) tname
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfend (idf)

      idv_curr = ncvid (idf, tname, iok)

      call ncvpt1 (idf, idv_curr, indx, tval, iok) 
      call ncsnc (idf, iok)

      idf_curr = idf

      return
      end

      subroutine odb_wrvar (idf, vname, var)
c17-------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname
      dimension ista(MAXVDIMS), icnt(MAXVDIMS), var(1)
      character*80 buff
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfend (idf)

      idv_curr = ncvid (idf, vname, iok) 
      call ncvinq (idf, idv_curr, buff, itype, ndim, icnt, nattr, iok) 

      do i = 1, ndim
         ista(i) = 1
         idd = icnt(i)
         call ncdinq (idf, idd, buff, icnt(i), iok)
      enddo

      call ncvpt (idf, idv_curr, ista, icnt, var, iok) 
      call ncsnc (idf, iok)

      idf_curr = idf

      return
      end

      subroutine odb_wrxv (idf, igrn, vname, var)
c18-------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname
      dimension ista(MAXVDIMS), icnt(MAXVDIMS), igrn(1), var(1)
      character*80 buff
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call odb_dfend (idf)

      idv_curr = ncvid (idf, vname, iok) 
      call ncvinq (idf, idv_curr, buff, itype, ndim, icnt, nattr, iok) 

      do i = 1, ndim
         if (igrn(i) .eq. 0) then
            ista(i) = 1
            idg = icnt(i)
            call ncdinq (idf, idg, buff, icnt(i), iok)
         else
            ista(i) = igrn(i)
            icnt(i) = 1
         endif
      enddo

      call ncvpt (idf, idv_curr, ista, icnt, var, iok) 
      call ncsnc (idf, iok)

      idf_curr = idf

      return
      end

      subroutine odb_wr1v2 (idf, igr2, vname, var)
c19-------------------------------------------------
      character*(*) vname
      dimension indx(2), var(1)

      indx(1) = 0
      indx(2) = igr2
      call odb_wrxv (idf, indx, vname, var)

      return
      end

      subroutine odb_wr1v3 (idf, igr2, igr3, vname, var)
c20--------------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = 0
      indx(2) = igr2
      indx(3) = igr3
      call odb_wrxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_wr1v4 (idf, igr2, igr3, igr4, vname, var)
c21--------------------------------------------------------------
      character*(*) vname
      dimension indx(4), var(1)

      indx(1) = 0
      indx(2) = igr2
      indx(3) = igr3
      indx(4) = igr4
      call odb_wrxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_wr2v3 (idf, igr3, vname, var)
c22--------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = 0
      indx(2) = 0
      indx(3) = igr3
      call odb_wrxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_wr2v4 (idf, igr3, igr4, vname, var)
c23-------------------------------------------------------
      character*(*) vname
      dimension indx(4), var(1)

      indx(1) = 0
      indx(2) = 0
      indx(3) = igr3
      indx(4) = igr4
      call odb_wrxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_wr3v4 (idf, igr4, vname, var)
c24--------------------------------------------------
      character*(*) vname
      dimension indx(4), var(1)

      indx(1) = 0
      indx(2) = 0
      indx(3) = 0
      indx(4) = igr4
      call odb_wrxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_rddm (idf, gn1, nn)
c25---------------------------------------
      character*(*) gn1
      character*80 buff
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      idd = ncdid (idf, gn1, iok)
      call ncdinq (idf, idd, buff, nn, iok)
      
      idf_curr = idf

      return
      end

      subroutine odb_rdgr (idf, gn1, nn, xx)
c26-------------------------------------------
      character*(*) gn1
      dimension xx(1), ista(1), icnt(1)
      character*80 buff
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      idd = ncdid (idf, gn1, iok)
      call ncdinq (idf, idd, buff, nn, iok)

      idv_curr = ncvid (idf, gn1, iok) 
      ista(1) = 1
      icnt(1) = nn
      call ncvgt (idf, idv_curr, ista, icnt, xx, iok) 

      idf_curr = idf

      return
      end

      subroutine odb_rdvar (idf, vname, var)
c27-------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname
      dimension ista(MAXVDIMS), icnt(MAXVDIMS), var(1)
      character*80 buff
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      idv_curr = ncvid (idf, vname, iok) 
      call ncvinq (idf, idv_curr, buff, itype, ndim, icnt, nattr, iok) 

      do i = 1, ndim
         ista(i) = 1
         idd = icnt(i)
         call ncdinq (idf, idd, buff, icnt(i), iok)
      enddo

      call ncvgt (idf, idv_curr, ista, icnt, var, iok) 

      idf_curr = idf

      return
      end

      subroutine odb_rdxv (idf, igrn, vname, var)
c28-------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      integer igrn(*)
      character*(*) vname
      dimension ista(MAXVDIMS), icnt(MAXVDIMS), var(1)
      character*80 buff
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      idv_curr = ncvid (idf, vname, iok) 
      call ncvinq (idf, idv_curr, buff, itype, ndim, icnt, nattr, iok) 

      do i = 1, ndim
         if (igrn(i) .eq. 0) then
            ista(i) = 1
            idg = icnt(i)
            call ncdinq (idf, idg, buff, icnt(i), iok)
         else
            ista(i) = igrn(i)
            icnt(i) = 1
         endif
      enddo

      call ncvgt (idf, idv_curr, ista, icnt, var, iok) 
      idf_curr = idf

      return
      end

      subroutine odb_rd1v2 (idf, igr2, vname, var)
c29-------------------------------------------------
      character*(*) vname
      dimension indx(2), var(1)

      indx(1) = 0
      indx(2) = igr2
      call odb_rdxv (idf, indx, vname, var)

      return
      end

      subroutine odb_rd1v3 (idf, igr2, igr3, vname, var)
c30--------------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = 0
      indx(2) = igr2
      indx(3) = igr3
      call odb_rdxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_rd1v4 (idf, igr2, igr3, igr4, vname, var)
c31--------------------------------------------------------------
      character*(*) vname
      dimension indx(4), var(1)

      indx(1) = 0
      indx(2) = igr2
      indx(3) = igr3
      indx(4) = igr4
      call odb_rdxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_rd2v3 (idf, igr3, vname, var)
c32--------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = 0
      indx(2) = 0
      indx(3) = igr3
      call odb_rdxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_rd2v4 (idf, igr3, igr4, vname, var)
c33-------------------------------------------------------
      character*(*) vname
      dimension indx(4), var(1)

      indx(1) = 0
      indx(2) = 0
      indx(3) = igr3
      indx(4) = igr4
      call odb_rdxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_rd3v4 (idf, igr4, vname, var)
c34--------------------------------------------------
      character*(*) vname
      dimension indx(4), var(1)

      indx(1) = 0
      indx(2) = 0
      indx(3) = 0
      indx(4) = igr4
      call odb_rdxv (idf, indx, vname, var)
      
      return
      end

      subroutine odb_getvdim (idf, vname, ndim, mdim)
c35-----------------------------------------------------
      character*(*) vname
      character*80 buff
      dimension mdim(1)

      idv_curr = ncvid (idf, vname, iok) 
      call ncvinq (idf, idv_curr, buff, itype, ndim, mdim, nattr, iok) 

      do i = 1, ndim
         idg = mdim(i)
         call ncdinq (idf, idg, buff, mdim(i), iok)
      enddo

      return
      end

      integer function idv_by_name (idf, vname)
c36--------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname
      if (vname .eq. 'global' .or. vname .eq. 'GLOBAL') then
         idv_by_name = NCGLOBAL
      else
         idv_by_name = ncvid (idf, vname, iok) 
      endif
      return
      end

      subroutine odb_setcattr (idf, vname, aname, val)
c371-------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname, val
      
      call odb_dfsta (idf)
      idv_curr = idv_by_name (idf, vname)
      call ncaptc (idf, idv_curr, aname, NCCHAR, lnblnk(val), val, iok)
      idf_curr = idf

      return
      end

      subroutine odb_getcattr (idf, vname, aname, val)
c372------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname, val
      
      call odb_dfsta (idf)
      idv_curr = idv_by_name (idf, vname)
      call ncainq (idf, idv_curr, aname, ttype, tlen, iok)
      if (tlen .gt. len(val)) then
         val = aname(1:len(aname))//' :too long!!!'
      else
         call ncagtc (idf, idv_curr, aname, val, tlen, iok)
      endif
      idf_curr = idf

      return
      end

      subroutine odb_setiattr (idf, vname, aname, ival)
c38-------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname

      call odb_dfsta (idf)
      idv_curr = idv_by_name (idf, vname)      
      call ncapt (idf, idv_curr, aname, NCLONG, 1, ival, iok)
      idf_curr = idf

      return
      end

      subroutine odb_setiaattr (idf, vname, aname, n, ival)
c39----------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname
      integer ival(1)

      call odb_dfsta (idf)
      idv_curr = idv_by_name (idf, vname)      
      call ncapt (idf, idv_curr, aname, NCLONG, n, ival, iok)
      idf_curr = idf

      return
      end

      subroutine odb_setrattr (idf, vname, aname, val)
c40-------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname

      call odb_dfsta (idf)
      idv_curr = idv_by_name (idf, vname)      
      call ncapt (idf, idv_curr, aname, NCFLOAT, 1, val, iok)
      idf_curr = idf

      return
      end

      subroutine odb_setraattr (idf, vname, aname, n, val)
c41---------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname
      real val(1)

      call odb_dfsta (idf)
      idv_curr = idv_by_name (idf, vname)      
      call ncapt (idf, idv_curr, aname, NCFLOAT, n, val, iok)
      idf_curr = idf

      return
      end

      subroutine odb_setdattr (idf, vname, aname, val)
c42-------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname
      double precision val

      call odb_dfsta (idf)
      idv_curr = idv_by_name (idf, vname)      
      call ncapt (idf, idv_curr, aname, NCDOUBLE, 1, val, iok)
      idf_curr = idf

      return
      end

      subroutine odb_getdattr (idf, vname, aname, val)
c43------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname
      double precision val

      idv_curr = idv_by_name (idf, vname)      
      call ncagt (idf, idv_curr, aname, val, iok)
      idf_curr = idf

      return
      end

      subroutine odb_getrattr (idf, vname, aname, val)
c44------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname

      idv_curr = idv_by_name (idf, vname)      
      call ncagt (idf, idv_curr, aname, val, iok)
      idf_curr = idf

      return
      end

      subroutine odb_getraattr (idf, vname, aname, n, val)
c45--------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname
      real val(1)

      idv_curr = idv_by_name (idf, vname)      
      call ncainq(idf, idv_curr, aname, itype, n, iok) 
      call ncagt (idf, idv_curr, aname, val, iok)
      idf_curr = idf

      return
      end

      subroutine odb_getiattr (idf, vname, aname, ival)
c46------------------------------------------------------
      character*(*) vname, aname      
      call odb_getrattr (idf, vname, aname, val)
      ival = int(val)
      return
      end

      subroutine odb_wrdesc (idf, desc)
c49-------------------------------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) desc

      call odb_dfsta (idf)
      call ncaptc (idf, NCGLOBAL, 'description', NCCHAR,lnblnk(desc),
     *             desc,iok) 

      return
      end

      subroutine odb_setfmt (fmt)
c50----------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr
      character*(*) fmt

      if (fmt(1:2).eq.'i1' .or. fmt(1:4).eq.'byte') then
         IDTYPE = NCBYTE
      elseif (fmt(1:4).eq.'char') then
         IDTYPE = NCCHAR
      elseif (fmt(1:2).eq.'i2' .or. fmt(1:5).eq.'short') then
         IDTYPE = NCSHORT
      elseif (fmt(1:2).eq.'i4' .or. fmt(1:4).eq.'long') then
         IDTYPE = NCLONG
      elseif (fmt(1:2).eq.'r4' .or. fmt(1:4).eq.'real' .or. 
     *        fmt(1:5).eq.'float') then
         IDTYPE = NCFLOAT
      elseif (fmt(1:2).eq.'r8' .or. fmt(1:6).eq.'double') then
         IDTYPE = NCDOUBLE
      else
         IDTYPE = NCFLOAT
      endif

      return
      end

      subroutine odb_setnxyt (idf,xname,nx,xx,yname,ny,yy,tname,nt,tt,
     *             nv,vnams)
c51-----------------------------------------------------------------------------
      character*(*) xname, yname, tname, vnams(1)
      dimension xx(1), yy(1), tt(1)
      logical odb_ifvar

      if (.not. odb_ifvar(idf, 'X') ) then 
         call odb_dfgr (idf, 'X', nx)
         call odb_wrgr (idf, 'X', xx)
         call odb_setcattr (idf, 'X', 'units', xname)
         call odb_setiattr (idf, 'X', 'gridtype', 0)
      endif

      if (.not. odb_ifvar(idf, 'Y') ) then 
         call odb_dfgr (idf, 'Y', ny)
         call odb_wrgr (idf, 'Y', yy)
         call odb_setcattr (idf, 'Y', 'units', yname)
         call odb_setiattr (idf, 'Y', 'gridtype', 0)
      endif

      if (.not. odb_ifvar(idf, 'T') ) then 
         call odb_dfgr (idf, 'T', nt)
         call odb_wrgr (idf, 'T', tt)
         call odb_setcattr (idf, 'T', 'units', tname)
         call odb_setiattr (idf, 'T', 'gridtype', 0)
      endif

      do iv = 1, nv
         call odb_dfvar3 (idf, 'X', 'Y', 'T', vnams(iv))
      enddo

      return
      end

      subroutine odb_setxyzt (idf,xname,nx,xx,yname,ny,
     *              yy,zname,nz,zz,tname,nt,tt,vnam)
c51-----------------------------------------------------------------------------
      character*(*) xname, yname, tname, zname, vnam
      dimension xx(1), yy(1), tt(1), zz(1)
      logical odb_ifvar

      if (.not. odb_ifvar(idf, 'X') ) then 
         call odb_dfgr (idf, 'X', nx)
         call odb_wrgr (idf, 'X', xx)
         call odb_setcattr (idf, 'X', 'units', xname)
         call odb_setiattr (idf, 'X', 'gridtype', 0)
      endif

      if (.not. odb_ifvar(idf, 'Y') ) then 
         call odb_dfgr (idf, 'Y', ny)
         call odb_wrgr (idf, 'Y', yy)
         call odb_setcattr (idf, 'Y', 'units', yname)
         call odb_setiattr (idf, 'Y', 'gridtype', 0)
      endif

      if (.not. odb_ifvar(idf, 'Z') ) then 
         call odb_dfgr (idf, 'Z', nz)
         call odb_wrgr (idf, 'Z', zz)
         call odb_setcattr (idf, 'Z', 'units', zname)
         call odb_setiattr (idf, 'Z', 'gridtype', 0)
      endif

      if (.not. odb_ifvar(idf, 'T') ) then 
         call odb_dfgr (idf, 'T', nt)
         call odb_wrgr (idf, 'T', tt)
         call odb_setcattr (idf, 'T', 'units', tname)
         call odb_setiattr (idf, 'T', 'gridtype', 0)
      endif

      call odb_dfvar4 (idf, 'X', 'Y', 'Z', 'T', vnam)

      return
      end

      subroutine odb_wrxy (idf, indxz, vname, var)
c62------------------------------------------------
      character*(*) vname
      dimension var(1)
      call odb_wr2v3 (idf, indxz, vname, var)
      return
      end

      subroutine odb_wrxz (idf, indxy, vname, var)
c63------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = 0
      indx(2) = indxy
      indx(3) = 0
      call odb_wrxv (idf, indx, vname, var)

      return
      end

      subroutine odb_wryz (idf, indxx, vname, var)
c64------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = indxx
      indx(2) = 0
      indx(3) = 0
      call odb_wrxv (idf, indx, vname, var)

      return
      end

      subroutine odb_wrxye (idf, time, vname, var)
c65------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = 0
      indx(2) = 0
      call odb_rddm (idf, 'T', ntdim)
      indx(3) = ntdim + 1

      call odb_wrxv (idf, indx, vname, var)
      call odb_wrtm (idf, 'T', indx(3), time)

      return
      end

      subroutine odb_wr2xye (idf, time, vn1, vn2, v1,v2)
c66------------------------------------------------------
      character*(*) vn1, vn2
      dimension indx(3), v1(1), v2(1)

      indx(1) = 0
      indx(2) = 0
      call odb_rddm (idf, 'T', ntdim)
      indx(3) = ntdim + 1

      call odb_wrxv (idf, indx, vn1, v1)
      call odb_wrxv (idf, indx, vn2, v2)
      call odb_wrtm (idf, 'T', indx(3), time)

      return
      end

      subroutine odb_rdxy (idf, indxz, vname, var)
c67------------------------------------------------
      character*(*) vname
      dimension var(1)
      call odb_rd2v3 (idf, indxz, vname, var)
      return
      end

      subroutine odb_rdxz (idf, indxy, vname, var)
c68------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = 0
      indx(2) = indxy
      indx(3) = 0
      call odb_rdxv (idf, indx, vname, var)

      return
      end

      subroutine odb_rdyz (idf, indxx, vname, var)
c69------------------------------------------------
      character*(*) vname
      dimension indx(3), var(1)

      indx(1) = indxx
      indx(2) = 0
      indx(3) = 0
      call odb_rdxv (idf, indx, vname, var)

      return
      end

      logical function odb_ifvar(idf, vname)
c70-------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      call ncpopt (0)
      idv_curr = ncvid (idf, vname, iok) 
      call ncpopt (NCVERBOS+NCFATAL)

      idf_curr = idf

      odb_ifvar = iok.eq.0 .and. idv_curr.ne.-1

      return
      end

      logical function odb_ifatt(idf, vname, aname)
c71--------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      character*(*) vname, aname
      common /odb_loc1/ IDTYPE, idf_curr, idv_curr

      idv_curr = idv_by_name (idf, vname)      
      call ncpopt (0)
      call ncainq (idf, idv_curr, aname, itype, ilen, iok)
      call ncpopt (NCVERBOS+NCFATAL)

      idf_curr = idf

      odb_ifatt = iok.eq.0

      return
      end

      subroutine odb_setxyt (idf,xname,nx,xx,yname,ny,yy,tname,
     *          nt,tt,vname)

      character*(*) xname, yname, tname, vname
      dimension xx(1), yy(1), tt(1)
      call odb_setnxyt (idf,xname,nx,xx,yname,ny,yy,tname,nt,tt,1,vname)
      return
      end


      subroutine odb_wipt (idf, ig, vname, var)
c----------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      integer idf, ivd, ierr, ig
      integer var
      character*(*) vname

      ivd = ncvid (idf, vname, ierr)
      call odb_setfmt('i4')
      call ncendf(idf,ierr)
      call ncvpt1 (idf, ivd, ig, var, ierr)
      call odb_setfmt('r4')
      call ncredf(idf,ierr)

      return
      end


      subroutine odb_wrpt (idf, ig, vname, var)
c----------------------------------------------------
#ifdef INC_PATH
      include 'INC_PATH'
#else
      include 'netcdf.inc'
#endif
      integer idf, ig(1), ivd, ierr
      real var
      character*(*) vname

      ivd = ncvid (idf, vname, ierr)
      call ncendf(idf,ierr)
      call ncvpt1 (idf, ivd, ig, var, ierr)
      call ncredf(idf,ierr)

      return
      end


