!  Wrapper for i/o buffer flushing routine
!  Author: George Seibel
!  Rewritten by: Meng-Juei Hsieh
!  Working for most Unix (BSD, Convex, Sun, Stellar, SGI Iris...)
subroutine amflsh(filenum)
   implicit none
   integer filenum ! unit file number to flush
   integer istat   ! return status from flush
#if defined(AIX) || defined( XLF90)
   call flush_(filenum) !page 222 in the V2.3 Language Reference of XLF
#else
#  ifdef SGI
      call flush(filenum,istat)
#  else
      call flush(filenum)
#  endif
#endif
   return
end subroutine amflsh