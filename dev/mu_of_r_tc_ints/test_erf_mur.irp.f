program test
 implicit none
 constant_mu = .True.
 touch constant_mu
! my_grid_becke = .False.
! touch my_grid_becke
 print*,'grid_type_sgn',grid_type_sgn
 call routine

end

subroutine routine
 implicit none
 double precision, allocatable :: ints_erf(:)
 integer :: i,j,k,l
 double precision, allocatable :: big_array_1(:,:,:,:),big_array_2(:,:,:,:)
 double precision :: exact, num_1, num_2,accu1, accu2, contrib1, contrib2 ,accu1_relat, accu2_relat
 double precision :: accu_ao
 allocate(ints_erf(ao_num), big_array_1(ao_num,ao_num,ao_num,ao_num),big_array_2(ao_num,ao_num,ao_num,ao_num) )
 big_array_1 = 0.d0
 call all_erf_mu_r1_lr_int_big_mat(big_array_1) ! erf(mu(r1)r12)/r12 mixed analytical/numerical
 big_array_2 = 0.d0
 call all_erf_mu_r1_lr_int_big_mat_bis(big_array_2)  ! exact 1/r12 + erfc(mu(r1)r12)/r12 mixed analytical/numerical
 accu_ao = 0.d0
 do l = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    call get_ao_two_e_integrals_erf(j,k,l,ao_num,ints_erf)
    do i = 1, ao_num
     exact = ints_erf(i)
     num_1 = big_array_1(i,k,j,l)
     num_2 = big_array_2(i,k,j,l)
     contrib1 = dabs(exact - num_1)
     contrib2 = dabs(exact - num_2)

     if(dabs(exact).lt.1.d-12)cycle
     accu_ao += 1.d0
     accu1_relat += contrib1/dabs(exact)
     accu1 += contrib1
     if(contrib1 .gt. 1.d-12)then
      print*,''
      print*,'num1 '
      print*,'i,k,j,l',i,k,j,l
      print*,'exact,num_1,difference'
      print*,exact,num_1,contrib1
      print*,''
     endif

     accu2_relat += contrib2/dabs(exact)
     accu2 += contrib2
     if(contrib2 .gt. 1.d-12)then
      print*,''
      print*,'num2 '
      print*,'i,k,j,l',i,k,j,l
      print*,'exact,num_2,difference'
      print*,exact,num_2,contrib2
      print*,''
     endif
     if(contrib2.gt.contrib1)then
      print*,''
      print*,'contrib2.gt.contrib1 !!'
      print*,'i,k,j,l',i,k,j,l
      print*,'contrib1,contrib2'
      print*,contrib1,contrib2
      print*,''
     endif
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,'accu1       = ',accu1/accu_ao
 print*,'accu1_relat = ',accu1_relat/accu_ao
 print*,''
 print*,'accu2       = ',accu2/accu_ao
 print*,'accu2_relat = ',accu2_relat/accu_ao



end
