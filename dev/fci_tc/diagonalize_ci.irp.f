
subroutine diagonalize_CI_dressed(ndet, E_tc,norm,pt2_data,print_pt2)
  use selection_types
  implicit none
  integer, intent(inout)              :: ndet      ! number of determinants from before 
  double precision, intent(inout)     :: E_tc,norm ! E and norm from previous wave function 
  type(pt2_type)  , intent(in)        :: pt2_data  ! PT2 from previous wave function 
  logical, intent(in) :: print_pt2
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
  print*,'*****'
  print*,'New wave function information'
  print*,'N_det tc               = ',N_det
  print*,'norm_ground_left_right = ',norm_ground_left_right
  print*,'eigval_right_tc = ',eigval_right_tc(1)
  print*,'Ndet, E_tc = ',N_det,eigval_right_tc(1)
  print*,'*****'
  if(print_pt2)then
   print*,'*****'
   print*,'previous wave function info'
   print*,'norm(before)      = ',norm
   print*,'E(before)         = ',E_tc
   print*,'PT1 norm          = ',dsqrt(pt2_data % overlap(1,1))
   print*,'E(before) + PT2   = ',E_tc + (pt2_data % pt2(1))/norm
   print*,'PT2               = ',pt2_data % pt2(1)
   print*,'Ndet, E_tc, E+PT2 = ',ndet,E_tc,E_tc + (pt2_data % pt2(1))/norm,dsqrt(pt2_data % overlap(1,1))
   print*,'*****'
  endif
  E_tc  = eigval_right_tc(1)
  if(cipsi_tc == "e_sym")then
   norm = norm_ground_right
  else
   norm  = norm_ground_left_right
  endif
  ndet  = N_det
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = reigvec_tc(i,j)
    enddo
  enddo
  do j=1,N_states
    do i=1,min(N_det,n_det_max_full)
      psi_left_guess(i,j) = leigvec_tc(i,j)
    enddo
  enddo
  SOFT_TOUCH psi_coef reigvec_tc leigvec_tc eigval_right_tc eigval_left_tc psi_left_guess
!  SOFT_TOUCH psi_coef 
  call routine_save_right 
end

subroutine print_CI_dressed(ndet, E_tc,norm,pt2_data,print_pt2)
  use selection_types
  implicit none
  integer, intent(inout)              :: ndet      ! number of determinants from before 
  double precision, intent(inout)     :: E_tc,norm ! E and norm from previous wave function 
  type(pt2_type)  , intent(in)        :: pt2_data  ! PT2 from previous wave function 
  logical, intent(in) :: print_pt2
  BEGIN_DOC
!  Replace the coefficients of the CI states by the coefficients of the
!  eigenstates of the CI matrix
  END_DOC
  integer :: i,j
  print*,'*****'
  print*,'New wave function information'
  print*,'N_det tc               = ',N_det
  print*,'norm_ground_left_right = ',norm_ground_left_right
  print*,'eigval_right_tc = ',eigval_right_tc(1)
  print*,'Ndet, E_tc = ',N_det,eigval_right_tc(1)
  print*,'*****'
  if(print_pt2)then
   print*,'*****'
   print*,'previous wave function info'
   print*,'norm(before)      = ',norm
   print*,'E(before)         = ',E_tc
   print*,'PT1 norm          = ',dsqrt(pt2_data % overlap(1,1))
   print*,'E(before) + PT2   = ',E_tc + (pt2_data % pt2(1))/norm
   print*,'PT2               = ',pt2_data % pt2(1)
   print*,'Ndet, E_tc, E+PT2 = ',ndet,E_tc,E_tc + (pt2_data % pt2(1))/norm,dsqrt(pt2_data % overlap(1,1))
   print*,'*****'
  endif
  E_tc  = eigval_right_tc(1)
  if(cipsi_tc == "e_sym")then
   norm = norm_ground_right
  else
   norm  = norm_ground_left_right
  endif
  ndet  = N_det
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = reigvec_tc(i,j)
    enddo
  enddo
  do j=1,N_states
    do i=1,min(N_det,n_det_max_full)
      psi_left_guess(i,j) = leigvec_tc(i,j)
    enddo
  enddo
  SOFT_TOUCH psi_coef reigvec_tc leigvec_tc eigval_right_tc eigval_left_tc psi_left_guess
end
