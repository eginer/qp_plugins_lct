program deriv_r12_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
 call test_all_prod_in_r
 call big_thing
 call big_thing_mo
!provide ao_two_e_eff_dr12_pot_array
!integer :: i,j,k,l
!double precision :: accu
!accu = 0.d0
!do j = 1, ao_num ! r2
! do i = 1, ao_num ! r1
!  do l = 1, ao_num ! r2 
!   do k = 1, ao_num ! r1 
!    accu += dabs(ao_two_e_eff_dr12_pot_array_no_cycle(k,l,i,j) - ao_two_e_eff_dr12_pot_array(k,l,i,j))
!   enddo 
!  enddo 
! enddo 
!enddo
!print*,'accu = ',accu/dble(ao_num)**4.d0
!
end
