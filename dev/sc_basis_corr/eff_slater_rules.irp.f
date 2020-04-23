
subroutine get_eff_mat_diag(det_i,h_eff_one_e,h_eff_two_e)
  BEGIN_DOC
! returns the diagonal matrix element of the effective hamiltonian coming from DFT 
  END_DOC
  use bitmasks
  implicit none
  integer(bit_kind), intent(in)  :: det_i(N_int,2)
  double precision, intent(out)  :: h_eff_one_e,h_eff_two_e
  integer                        :: n_occ_ab(2)
  integer                        :: occ(N_int*bit_kind_size,2)
  integer :: i,j,h1,h2,p1,p2
  double precision :: eff_int_mat(mo_num,mo_num)
  double precision               :: get_two_e_integral,integral,integral_bis,get_mo_two_e_int_mu_of_r
  provide tot_eff_two_e
  call bitstring_to_list_ab(det_i, occ, n_occ_ab, N_int)
  h_eff_one_e = 0.d0
  h_eff_two_e = 0.d0
  ! alpha <-> alpha 
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   h_eff_one_e += 0.5d0 * ( pot_basis_alpha_mo_su_pbe_ot(h1,h1,1) + pot_basis_beta_mo_su_pbe_ot(h1,h1,1) )
  enddo

  ! beta  <-> beta  
  do i = 1, n_occ_ab(2)
   h1 = occ(i,2)
   h_eff_one_e += 0.5d0 * ( pot_basis_alpha_mo_su_pbe_ot(h1,h1,1) + pot_basis_beta_mo_su_pbe_ot(h1,h1,1) )
  enddo

  ! alpha <-> beta  
  do i = 1, n_occ_ab(1)
   h1 = occ(i,1)
   do j = 1, n_occ_ab(2)
    h2 = occ(j,2)
    integral = eff_two_e(h1,h2,h1,h2,1)
    h_eff_two_e += integral
   enddo
  enddo

end

subroutine get_eff_mat_off_diag(det_i,det_j,h_eff_one_e,h_eff_two_e)
  BEGIN_DOC
! returns the OFF-diagonal matrix element of the effective hamiltonian coming from DFT 
  END_DOC
  use bitmasks
  implicit none
  integer(bit_kind), intent(in)  :: det_i(N_int,2),det_j(N_int,2)
  double precision, intent(out)  :: h_eff_one_e,h_eff_two_e
  integer                        :: n_occ_ab(2),exc(0:2,2,2)
  integer                        :: occ(N_int*bit_kind_size,2)
  integer :: i,j,h1,h2,p1,p2,degree,spin,other_spin(2)
  double precision               :: get_two_e_integral,integral,integral_bis
  double precision               :: hij,get_mo_two_e_int_mu_of_r,phase
  other_spin(1) = 2
  other_spin(2) = 1
  h_eff_two_e = 0.d0
  h_eff_one_e = 0.d0
  call get_excitation_degree(det_i,det_j,degree,N_int)
  if(degree == 2)then
   call get_double_excitation(det_i,det_j,exc,phase,N_int)
   if (exc(0,1,1) == 1) then
    ! Single alpha, single beta
    h1 = exc(1,1,1)
    h2 = exc(1,1,2)
    p1 = exc(1,2,1)
    p2 = exc(1,2,2) 
    integral = eff_two_e(h1,h2,p1,p2,1)
    h_eff_two_e = phase * integral
   endif
  else 
   call get_single_excitation(det_i,det_j,exc,phase,N_int)
   call bitstring_to_list_ab(det_i, occ, n_occ_ab, N_int)
   if (exc(0,1,1) == 1) then
    ! Single alpha
    h1 = exc(1,1,1)
    p1 = exc(1,2,1)
    spin = 1
   else
    ! Single beta
    h1 = exc(1,1,2)
    p1 = exc(1,2,2)
    spin = 2
   endif
   h_eff_one_e = 0.5d0 * phase * ( pot_basis_alpha_mo(h1,p1,1) + pot_basis_beta_mo(h1,p1,1) )
   do i = 1, n_occ_ab(other_spin(spin))
    h2 = occ(i,other_spin(spin))
    p2 = h2
    integral = eff_two_e(h1,h2,p1,p2,1)
    h_eff_two_e += integral * phase
   enddo
  endif
end


subroutine i_eff_H_j(det_i,det_j,h_eff_one_e,h_eff_two_e)
  BEGIN_DOC
! returns the general matrix element of the effective hamiltonian coming from DFT 
  END_DOC
  use bitmasks
  implicit none
  integer(bit_kind), intent(in)  :: det_i(N_int,2),det_j(N_int,2)
  double precision, intent(out)  :: h_eff_one_e,h_eff_two_e
  integer :: degree
  h_eff_one_e = 0.d0
  h_eff_two_e = 0.d0
  call get_excitation_degree(det_i,det_j,degree,N_int)
  if(degree.gt.2)then
   return 
  else if(degree == 0)then
   call get_eff_mat_diag(det_i,h_eff_one_e,h_eff_two_e)
  else 
   call get_eff_mat_off_diag(det_i,det_j,h_eff_one_e,h_eff_two_e)
  endif
 end
