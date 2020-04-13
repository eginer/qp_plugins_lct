program test
 implicit none
 call routine

end

subroutine routine
 implicit none
 integer :: i,j,k,l
 double precision :: mo_integrals(mo_num,mo_num)
 double precision :: accu
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   call compute_all_ijkl_for_jl_mu_of_r_int(j,i,mo_integrals)
   do k = 1, mo_num
    do l = 1, mo_num
     if(dabs(eff_two_e(l,k,j,i,1)).gt.1.d-10)then
      if(dabs(eff_two_e(l,k,j,i,1) - mo_integrals(l,k)).gt.1.d-10)then
       print*,eff_two_e(l,k,j,i,1), mo_integrals(l,k),dabs(eff_two_e(l,k,j,i,1) - mo_integrals(l,k))
      endif
     endif
     accu += dabs(eff_two_e(l,k,j,i,1) - mo_integrals(l,k))
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num**4)
end
