subroutine test_new_erf_ints
 implicit none
 integer :: i,j,k,l
 double precision :: ao_two_e_integral_schwartz_accel_erf_new,ao_two_e_integral_schwartz_accel_erf
 double precision :: ref,new
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = 1, ao_num
     ref = ao_two_e_integral_schwartz_accel_erf(i,j,k,l)
     new = ao_two_e_integral_schwartz_accel_erf_new(i,j,k,l)
     if(dabs(ref - new).gt.1.d-14)then
       print*,'Problemmmmmm'
       print*,i,j,k,l
       stop
     endif
    enddo
   enddo
  enddo
 enddo

end
