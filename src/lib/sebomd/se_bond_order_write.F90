subroutine se_bond_order_write(natoms, iatnum, natorb, ip1, ijmat, pdiat)
  use file_io_dat, only : sebounit
  use sebomd_module, only : sebomd_obj
  implicit none
  integer :: natoms
  integer, dimension (*) :: iatnum, natorb, ip1, ijmat
  double precision, dimension (*) :: pdiat

  integer :: i,j
  integer :: iat,jat
  integer :: norbsi,norbsj
  integer :: npairs
  integer :: ij, ijaddr, ijstart, ijend
  double precision :: bk
  integer :: k
  double precision :: bocut

  bocut = sebomd_obj%bocut
  write(sebounit,'("BOND ORDER natoms =",i6, " bocut = ",f8.5)') natoms, bocut
  npairs = ip1(natoms+1)-1
  do i = 2, natoms
    iat = iatnum(i)
    if (iat.ne.0) then
      norbsi = natorb(iat)
      do j = 1, i-1
        jat = iatnum(j)
        if (jat.ne.0) then
          norbsj = natorb(jat)
          call se_ijfind(npairs,i,j,ijaddr)
          if (ijaddr.ne.0) then
            ijstart = ijmat(ijaddr)
            ijend = ijstart + norbsi*norbsj -1
            bk = 0.0d0
            do ij=ijstart,ijend
              bk = bk + pdiat(ij)*pdiat(ij)
            enddo
            if (bk.gt.bocut) then
              write(sebounit,'(2i5,f8.3)')  j,i,bk
            endif
          endif
        endif
      end do
    endif
  end do
end subroutine se_bond_order_write
