
! ---

BEGIN_PROVIDER [ double precision, good_hermit_tc_fock_mat, (mo_num, mo_num)]

  implicit none
  integer :: i, j

  good_hermit_tc_fock_mat = Fock_matrix_tc_mo_tot
  do j = 1, mo_num
    do i = 1, j-1
      good_hermit_tc_fock_mat(i,j) = Fock_matrix_tc_mo_tot(j,i) 
    enddo
  enddo

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, grad_good_hermit_tc_fock_mat]

  implicit none
  integer :: i, j

  grad_good_hermit_tc_fock_mat = 0.d0
  do i = 1, elec_alpha_num
    do j = elec_alpha_num+1, mo_num
      grad_good_hermit_tc_fock_mat += dabs(good_hermit_tc_fock_mat(i,j))
    enddo
  enddo

END_PROVIDER 

! ---

subroutine save_good_hermit_tc_eigvectors()

  implicit none
  integer        :: sign
  character*(64) :: label
  logical        :: output

  sign = 1
  label = "Canonical"
  output = .False.
  
  call mo_as_eigvectors_of_mo_matrix(good_hermit_tc_fock_mat, mo_num, mo_num, label, sign, output)                                                    
end subroutine save_good_hermit_tc_eigvectors

! ---

