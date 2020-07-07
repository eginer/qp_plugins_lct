BEGIN_PROVIDER [double precision, ao_deriv_r12_ints, (ao_num,ao_num,ao_num,ao_num)]
 implicit none
 integer :: i,j,k,l
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num

    enddo
   enddo
  enddo
 enddo
END_PROVIDER 
