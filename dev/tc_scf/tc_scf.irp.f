program tc_scf

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call routine_scf()

end

subroutine routine_scf()

  implicit none
  integer          :: i, j, it
  double precision :: e_save, e_delta

  it = 0
  print*,'iteration = ', it

  !print*,'grad_good_hermit_tc_fock_mat = ', grad_good_hermit_tc_fock_mat
  print*,'***'
  print*,'TC HF total energy = ', TC_HF_energy
  print*,'TC HF 1 e   energy = ', TC_HF_one_electron_energy
  print*,'TC HF 2 e   energy = ', TC_HF_two_e_energy
  if(.not. bi_ortho)then
  print*,'TC HF 3 body       = ', diag_three_elem_hf
  endif
  print*,'***'

  e_delta = 10.d0
  e_save  = TC_HF_energy

  if(bi_ortho)then
   mo_l_coef = fock_tc_leigvec_ao
   mo_r_coef = fock_tc_reigvec_ao
   call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
   call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
   TOUCH mo_l_coef mo_r_coef
  else
   print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
   call save_good_hermit_tc_eigvectors
   TOUCH mo_coef 
  endif

  if(bi_ortho)then
   do while( it .lt. n_it_tcscf_max )
     it += 1
     print*,'iteration = ', it
     print*,'***'
     print*,'TC HF total energy = ', TC_HF_energy
     print*,'TC HF 1 e   energy = ', TC_HF_one_electron_energy
     print*,'TC HF 2 non hermit = ', TC_HF_two_e_energy
!     print*,'TC HF 3 body       = ', diag_three_elem_hf
     print*,'***'
     e_delta = dabs( TC_HF_energy - e_save )
     print*, 'it, delta E = ', it, e_delta
     e_save = TC_HF_energy
     mo_l_coef = fock_tc_leigvec_ao
     mo_r_coef = fock_tc_reigvec_ao
     call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
     call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
     TOUCH mo_l_coef mo_r_coef
   enddo
  else
   do while( (grad_good_hermit_tc_fock_mat.gt.thresh_tcscf) .or. it .lt. n_it_tcscf_max )
      print*,'grad_good_hermit_tc_fock_mat = ',grad_good_hermit_tc_fock_mat
      it += 1
      print*,'iteration = ', it
      print*,'***'
      print*,'TC HF total energy = ', TC_HF_energy
      print*,'TC HF 1 e   energy = ', TC_HF_one_electron_energy
      print*,'TC HF 2 e   energy = ', TC_HF_two_e_energy
      print*,'TC HF 3 body       = ', diag_three_elem_hf
      print*,'***'
      call save_good_hermit_tc_eigvectors
      TOUCH mo_coef 
   enddo
  endif

end subroutine routine_scf

! ---

