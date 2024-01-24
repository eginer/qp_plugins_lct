program print_tc_lr_mos
 implicit none
 integer :: i,j
 print*,'left-coef'
 do i = 1, ao_num
  write(*,'(1000(F10.6,X))')mo_l_coef(i,:)
 enddo

end
