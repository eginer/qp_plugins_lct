
subroutine diagonalize_CI_dressed
  implicit none
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = reigvec_tc(i,j)
    enddo
  enddo
  SOFT_TOUCH psi_coef
end
