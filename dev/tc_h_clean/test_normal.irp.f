program test
 implicit none
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call test_pouet
end

subroutine test_pouet
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), allocatable :: det_i(:,:)
 allocate(det_i(N_int,2))
 integer :: h1,p1,s1,s2,i_ok,i,h2,p2
 double precision :: accu(2,2)
 integer   :: exc(0:2,2,2)
 integer   :: degree
 double precision :: phase,hmono,heff,hderiv,hthree,htot,hthree_new
 double precision :: hnew,accu_n(2,2)
 accu = 0.d0
 accu_n = 0.d0
! do h1 = 1, mo_num
!  do h2 = 1, mo_num
!   do p1 = h1, mo_num
!    do p2 = h2, mo_num
 do h1 = 1, 1
  do p1 = 7,7
   do h2 = 2,2 
    do p2 = 4,4
!     do s1 = 1, 2
!      do s2 = 1, 2
     do s1 = 1,1
      do s2 = 1,1
       do i = 1, N_int
         det_i(i,1) = ref_bitmask(i,1)
         det_i(i,2) = ref_bitmask(i,2)
       enddo
       call do_single_excitation(det_i,h1,p1,s1,i_ok)
       if(i_ok == -1)cycle
       call do_single_excitation(det_i,h2,p2,s2,i_ok)
       if(i_ok == -1)cycle
       call get_excitation_degree(ref_bitmask,det_i,degree,N_int)
       if(degree.ne.2)cycle
       accu_n(s2,s1) += 1.d0
       call get_excitation(ref_bitmask,det_i,exc,degree,phase,N_int)
       call htilde_mu_mat(ref_bitmask,det_i,hmono,heff,hderiv,hthree,htot)
       hthree_new = phase*(normal_two_body_aa_bb(p2,h2,p1,h1))
       hnew = hmono+heff+hderiv+hthree_new

!       if(dabs(hthree).lt.1.d-10)cycle
       if(dabs(htot-hnew).gt.1.d-6)then
!        print*,s1,s2,htot,hnew,dabs(hnew - htot)
        print*,h1,h2,p1,p2
        print*,s1,s2,hthree,hthree_new,dabs(hthree_new - hthree)
       endif
       accu(s2,s1) += dabs(hnew - htot)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,''
 do i = 1, 2
  print*,accu(1,i)/accu_n(1,i),accu(2,i)/accu_n(2,i)
 enddo
  
end
