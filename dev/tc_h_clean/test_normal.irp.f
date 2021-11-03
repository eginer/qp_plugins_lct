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
 double precision, allocatable :: array(:,:,:,:)
 allocate(array(mo_num,mo_num,mo_num,mo_num))
 array = 0.d0
 accu = 0.d0
 accu_n = 0.d0
 integer, allocatable           :: occ(:,:),Ne(:)
 allocate(occ(N_int*bit_kind_size,2),Ne(2))
 call bitstring_to_list_ab(ref_bitmask,occ,Ne,N_int)
! do h1 = 1, 1
!  do h2 = 2,2 
!   do p1 = 7,7
!    do p2 = 3,3
 do h1 = 1, mo_num
  do h2 = 1, mo_num
   do p1 = 1, mo_num
    do p2 = 1, mo_num
     do s1 = 1, 2
      do s2 = 1, 2
!     do s1 = 2,2
!      do s2 = 2,2
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
       call get_excitation(det_i,ref_bitmask,exc,degree,phase,N_int)
       call htilde_mu_mat(det_i,ref_bitmask,hmono,heff,hderiv,hthree,htot)
       integer :: hh1,pp1,hh2,pp2,ss1,ss2
       call decode_exc(exc,2,hh1,pp1,hh2,pp2,ss1,ss2)
!       call give_aab_contraction(hh1,hh2,pp1,pp2,Ne,occ,hthree_new)
!       print*,'phase = ',phase
!       print*,'array(p2,h2,p1,h1)',array(pp2,hh2,pp1,hh1)
       if(s1==s2)then
        hthree_new = normal_two_body_aa_bb(pp2,hh2,pp1,hh1)
       else
        hthree_new = normal_two_body_ab(pp2,hh2,pp1,hh1)
       endif
       hthree_new *= phase
       hnew = hmono+heff+hderiv+hthree_new
!
!       print*,hthree,hthree_new,dabs(hthree_new - hthree)
!       if(dabs(hthree).lt.1.d-10)cycle
       if(dabs(htot-hnew).gt.1.d-10)then
        print*,'***'
        print*,'wrong'
!        print*,s1,s2,htot,hnew,dabs(hnew - htot)
        print*,h1,h2,p1,p2,s1,s2
        print*,hthree,hthree_new,dabs(hthree_new - hthree)
        print*,'***'
        stop
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
