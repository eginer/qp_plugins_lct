
subroutine diagonalize_CI_dressed(E_tc,pt2_data,print_pt2)
  use selection_types
  implicit none
  double precision, intent(out)  :: E_tc
  type(pt2_type)  , intent(in)   :: pt2_data
  logical, intent(in) :: print_pt2
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
  print*,'N_det tc        = ',N_det
  print*,'eigval_right_tc = ',eigval_right_tc(1)
  if(print_pt2)then
   print*,'E+PT2           = ',eigval_right_tc(1) + pt2_data % pt2(1)
   print*,'PT2             = ',pt2_data % pt2(1)
  endif
  E_tc  = eigval_right_tc(1)
end
