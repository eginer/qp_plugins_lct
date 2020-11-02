program print_hcore
 implicit none
 integer :: i
 open(1, file = 'hcore') 
 do i = 1, mo_num
  write(1,'(100(F16.10,X))')mo_one_e_integrals(i,:)
 enddo
 close(1)
end
