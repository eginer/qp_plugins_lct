program tc_scf
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
  call routine_mo
end

subroutine routine
 implicit none
 integer :: i,j
 double precision :: accu_alpha,accu_beta
 accu_alpha = 0.D0
 accu_beta = 0.D0
 do i = 1, ao_num
  do j = 1, ao_num
   accu_alpha += dabs(two_e_tc_hermit_integral_alpha(j,i) - ao_two_e_integral_alpha(j,i))
   accu_beta  += dabs(two_e_tc_hermit_integral_beta(j,i) - ao_two_e_integral_beta(j,i))
  enddo
 enddo
 print*,'accu_beta  = ',accu_beta
 print*,'accu_alpha = ',accu_alpha

end

subroutine routine_mo
 implicit none
 integer :: i,j,i_ok
 double precision :: f_tc
 double precision :: accu_alpha,accu_beta,hmono,heff,hderiv,hthree,htot
 integer(bit_kind), allocatable :: det_i(:,:)
 allocate(det_i(N_int,2))
 do i = 1, elec_alpha_num
  do j = elec_alpha_num+1, mo_num
   det_i(:,1) = ref_bitmask(:,1)
   det_i(:,2) = ref_bitmask(:,2)
   f_tc = Fock_matrix_tc_mo_alpha(i,j) 
   call do_single_excitation(det_i,i,j,1,i_ok)
   call htilde_mu_mat(ref_bitmask,det_i,hmono,heff,hderiv,hthree,htot)
   print*,'i,j',i,j
   print*,'ref,new,dabs'
   print*,htot,f_tc, dabs(f_tc - htot)
   accu_alpha += (f_tc - htot)
  enddo
 enddo


end
