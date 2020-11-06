program print_hcore
 implicit none
 call test_one_e 

end

subroutine test_one_e
 implicit none
 integer :: i,j
 double precision :: accu
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   accu+= (one_e_dm_mo_alpha_for_dft(j,i,1) + one_e_dm_mo_alpha_for_dft(j,i,1)) * mo_one_e_integrals(j,i)
  enddo
 enddo
 print*,'accu = ',accu
end
