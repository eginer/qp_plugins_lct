program eff_two_e
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  double precision :: get_two_e_integral,get_mo_two_e_int_mu_of_r
  double precision :: integral,integral_2
  integer :: i,j,k,l
  i = 1
  j = 1
  k = 1
  l = 1
  double precision :: accu,tmp
  accu = 0.d0
  do i = 1, mo_num
   do j = 1, mo_num
    do k = 1, mo_num
     do l = 1, mo_num
      integral_2 = get_mo_two_e_int_mu_of_r(i,j,k,l,mo_int_mu_of_r_map)
      integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
      tmp = dabs(integral - integral_2)
      accu += tmp
!      if(dabs(integral).lt.1.d-10)cycle
!      if(tmp.gt.1.d-10)then
!       print*,'ahahahaha'
!       print*,i,j,k,l
!       print*,integral,integral_2,tmp
!      endif
     enddo
    enddo
   enddo
  enddo
  print*,'accu     ',accu
  print*,'accuav = ',accu/dble(mo_num**4)
end
