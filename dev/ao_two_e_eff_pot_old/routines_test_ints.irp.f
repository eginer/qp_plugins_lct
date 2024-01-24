subroutine test_new_erf_ints
 implicit none
 integer :: i,j,k,l
 double precision :: ao_two_e_integral_schwartz_accel_erf_new,ao_two_e_integral_schwartz_accel_erf
 double precision :: ref,new
 double precision :: wall0,wall1
 call wall_time(wall0)
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = j,j
     ref = ao_two_e_integral_schwartz_accel_erf(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'wall time old = ',wall1 - wall0
 double precision :: accu
 accu = 0.d0
 call wall_time(wall0)
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    do l = j,j
     ref = ao_two_e_integral_schwartz_accel_erf(i,j,k,l)
     new = ao_two_e_integral_schwartz_accel_erf_new(i,j,k,l)
     accu += dabs(ref - new)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'wall time new = ',wall1 - wall0
 print*,'accu = ',accu

end
