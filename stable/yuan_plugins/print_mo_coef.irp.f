program print_mos
 integer :: i
 open(1, file = 'mo_coef')                                                                                                               
 do i = 1, mo_num
  write(1,'(100(F16.10,X))')mo_coef(i,1:mo_num)
 enddo
 close(1) 


end
