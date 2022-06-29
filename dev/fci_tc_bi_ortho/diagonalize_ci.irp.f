
 BEGIN_PROVIDER [ double precision, psi_selectors_rcoef_bi_orth_transp, (N_states, psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_selectors_lcoef_bi_orth_transp, (N_states, psi_det_size) ]

  implicit none
  integer :: i, k

  do i = 1, N_det_selectors
    do k = 1, N_states
      psi_selectors_rcoef_bi_orth_transp(k,i) = reigvec_tc_bi_orth_sorted(i,k)
      psi_selectors_lcoef_bi_orth_transp(k,i) = leigvec_tc_bi_orth_sorted(i,k)
    enddo
  enddo

END_PROVIDER

subroutine diagonalize_CI_tc_bi_ortho(ndet, E_tc,norm,pt2_data,print_pt2)
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
  print*,'norm_ground_left_right_bi_orth = ',norm_ground_left_right_bi_orth
  print*,'eigval_right_tc = ',eigval_right_tc_bi_orth(1)
  print*,'Ndet, E_tc = ',N_det,eigval_right_tc_bi_orth(1)
  print*,'*****'
  if(print_pt2)then
   print*,'*****'
   print*,'previous wave function info'
   print*,'norm(before)      = ',norm
   print*,'E(before)         = ',E_tc + nuclear_repulsion
   print*,'PT1 norm          = ',dsqrt(pt2_data % overlap(1,1))
   print*,'E(before) + PT2   = ',E_tc + nuclear_repulsion + (pt2_data % pt2(1))/norm
   print*,'PT2               = ',pt2_data % pt2(1)
   print*,'Ndet, E_tc, E+PT2 = ',ndet,E_tc + nuclear_repulsion,E_tc + nuclear_repulsion + (pt2_data % pt2(1))/norm,dsqrt(pt2_data % overlap(1,1))
   print*,'*****'
  endif
  E_tc  = eigval_right_tc_bi_orth(1)
  norm  = norm_ground_left_right_bi_orth
  ndet  = N_det
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = reigvec_tc_bi_orth(i,j)
    enddo
  enddo
!  do j=1,N_states
!    do i=1,min(N_det,n_det_max_full)
!      psi_left_guess(i,j) = leigvec_tc_bi_orth(i,j)
!    enddo
!  enddo
  SOFT_TOUCH  eigval_left_tc_bi_orth  eigval_right_tc_bi_orth  leigvec_tc_bi_orth  norm_ground_left_right_bi_orth psi_coef reigvec_tc_bi_orth 

!  call routine_save_right 
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
  print*,'norm_ground_left_right_bi_orth = ',norm_ground_left_right_bi_orth
  print*,'eigval_right_tc = ',eigval_right_tc_bi_orth(1)
  print*,'Ndet, E_tc = ',N_det,eigval_right_tc_bi_orth(1)
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
  E_tc  = eigval_right_tc_bi_orth(1)
  norm  = norm_ground_left_right_bi_orth
  ndet  = N_det
  do j=1,N_states
    do i=1,N_det
      psi_coef(i,j) = reigvec_tc_bi_orth(i,j)
    enddo
  enddo
!  do j=1,N_states
!    do i=1,min(N_det,n_det_max_full)
!      psi_left_guess(i,j) = leigvec_tc_bi_orth(i,j)
!    enddo
!  enddo
  SOFT_TOUCH  eigval_left_tc_bi_orth  eigval_right_tc_bi_orth  leigvec_tc_bi_orth  norm_ground_left_right_bi_orth  psi_coef  reigvec_tc_bi_orth 

  !psi_left_guess
end
