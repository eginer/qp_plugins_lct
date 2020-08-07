program plugins_test
  implicit none
 integer :: i,j
 do i = 1, mo_num
  do j = 1, mo_num
!   write(34,*)i,j,potential_c_alpha_mo_basis_pbe_ueg(j,i,1)
  enddo
 enddo

end
