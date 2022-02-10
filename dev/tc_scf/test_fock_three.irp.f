program test_fock_three
 implicit none
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
! call routine_test
! call routine
! call test_direct
! call test_diag_three
 call test_scaled_fock_three

end

subroutine test_scaled_fock_three
 implicit none
 integer :: i,j
 double precision :: accu, ref, new,new_2
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   call give_contrib_three_fock(i,j,ref)
   call give_fock_ia_same_spin(i,j,new)
   call give_fock_ia_scaled_op_spin(i,j,new_2)
   new += new_2
   accu += dabs(ref - new)
   if(dabs(ref) .gt.1.d-10 )then
    print*,ref,new,ref/new
   endif
  enddo
 enddo
 print*,'accu = ',accu

end

subroutine test_diag_three
 implicit none
 double precision :: hmono,heff,hderiv,hthree,htot
 call htilde_mu_mat(ref_bitmask,ref_bitmask,hmono,heff,hderiv,hthree,htot)
 print*,' hthree                 ',hthree
 print*,' diag_three_elem_hf     ',diag_three_elem_hf
 print*,' diag_three_elem_hf_old ',diag_three_elem_hf_old
 print*,'****'
 print*,'htot                    ',htot
 print*,'TC_right_HF_energy      ',TC_right_HF_energy

end

subroutine test_direct
 implicit none
 integer :: i,a,k,l
 double precision :: ref, new,accu

 accu = 0.d0
 do i = 1, mo_num
  do a = 1, mo_num
   call give_contrib_three_fock(i,a,ref)
!   call give_fock_ia_real_space_old(i,a,ref)
!   call give_fock_ia_real_space_old_bis(i,a,new)
!   call give_fock_ia_real_space_new(i,a,new)
    call give_fock_ia_real_space_prov(i,a,new)
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
 print*,''
  print*,'***'
   print*,'TC HF total energy = ',TC_right_HF_energy
   print*,'TC HF 1 e   energy = ',TC_right_HF_one_electron_energy
   print*,'TC HF 2 e hermit   = ',TC_right_HF_two_e_hermit_energy
   print*,'TC HF 2 non hermit = ',TC_right_HF_two_e_n_hermit_energy
  print*,'***'
 print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
 allocate(det_i(N_int,2))
 do i = 1, elec_alpha_num
  do a = elec_alpha_num+1, mo_num
   det_i(:,1) = ref_bitmask(:,1)
   det_i(:,2) = ref_bitmask(:,2)
   call do_single_excitation(det_i,i,a,1,i_ok)
   call give_contrib_three_fock(i,a,f_tc)
!   f_tc = fock_3_mat(i,a) ! <HF|(a^dagger_a a_i) H |HF > = F(i,a)
   call htilde_mu_mat(det_i,ref_bitmask,hmono,heff,hderiv,hthree,htot)
!   print*,'***'
!   print*,'hthree = ',hthree
!   print*,'f_tc   = ',f_tc 
   if(dabs(hthree).gt.1.d-10)then
    print*,'i,a',i,a
    print*,'ref,new,dabs'
    print*,dabs(hthree),dabs(f_tc), dabs(f_tc/hthree)
   endif
   accu_alpha += dabs(dabs(f_tc) - dabs(hthree))
  enddo
 enddo
 print*,'accu_alpha = ',accu_alpha
end
