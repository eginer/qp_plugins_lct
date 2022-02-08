program test_fock_three
 implicit none
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
! call routine_test
! call routine
 call test_direct

end

subroutine test_direct
 implicit none
 integer :: i,a,k,l
 double precision :: ref, new,accu

 accu = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
!   call give_contrib_three_fock(i,a,ref)
   call give_fock_ia_real_space(i,a,ref)
   call give_fock_ia_real_space_new(i,a,new)
   accu += dabs(ref - new)
   if(dabs(ref).gt.1.d-10)then
    print*,'***'
    print*,i,a
    print*,ref,new, dabs(ref - new)
   endif
  enddo
 enddo
 print*,''
 print*,''
 print*,'accu = ',accu

end

subroutine routine_test
 implicit none
 integer :: i,a
 double precision :: accu
 accu = 0.d0
 do i = 1, elec_alpha_num
  do a = elec_alpha_num+1, mo_num
   accu += dabs(ref_fock_three_new(i,a) - ref_fock_three(i,a))
   if(dabs(ref_fock_three(i,a)).gt.1.d-10)then
    print*,'i,a',i,a
    print*,'ref_fock_three_new , ref_fock_three',ref_fock_three(i,a),ref_fock_three_new(i,a)
   endif
  enddo
 enddo
 print*,'accu = ',accu


end

subroutine routine
 implicit none
 integer :: i,a,i_ok
 double precision :: f_tc
 double precision :: accu_alpha,accu_beta,hmono,heff,hderiv,hthree,htot
 integer(bit_kind), allocatable :: det_i(:,:)
 print*,'e_tilde_00 = ',e_tilde_00
 print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
 allocate(det_i(N_int,2))
 do i = 1, elec_alpha_num
  do a = elec_alpha_num+1, mo_num
   det_i(:,1) = ref_bitmask(:,1)
   det_i(:,2) = ref_bitmask(:,2)
   call do_single_excitation(det_i,i,a,1,i_ok)
   f_tc = ref_fock_three_new(i,a) ! <HF|H a^dagger_a a_i |HF > = F(i,a)
   call htilde_mu_mat(det_i,ref_bitmask,hmono,heff,hderiv,hthree,htot)
   if(dabs(hthree).gt.1.d-10)then
    print*,'i,a',i,a
    print*,'ref,new,dabs'
    print*,dabs(hthree),dabs(f_tc), dabs(f_tc) - dabs(hthree)
   endif
   accu_alpha += dabs(dabs(f_tc) - dabs(hthree))
  enddo
 enddo
 print*,'accu_alpha = ',accu_alpha
end
