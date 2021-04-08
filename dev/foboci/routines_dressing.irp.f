subroutine provide_matrix_dressing(dressing_matrix,ndet_generators_input,psi_det_generators_input)
 use bitmasks
 implicit none
 integer, intent(in) :: ndet_generators_input
 integer(bit_kind), intent(in) :: psi_det_generators_input(N_int,2,ndet_generators_input)
 double precision, intent(inout) :: dressing_matrix(ndet_generators_input,ndet_generators_input)
 double precision :: H_array(N_det),hka
 logical :: is_a_ref_det(N_det)
 integer :: i,j,n_det_ref_tmp
 integer :: connected_to_ref_by_mono,degree
 double precision :: coef_ref(Ndet_generators_input)
 double precision :: accu,lambda_i
 integer :: k
 integer :: index_ref_tmp(N_det)
 is_a_ref_det = .False.
 n_det_ref_tmp = 0
 do i = 1, N_det
  do j = 1, Ndet_generators_input
   call get_excitation_degree(psi_det(1,1,i),psi_det_generators_input(1,1,j),degree,N_int)  
   if(degree == 0)then
    is_a_ref_det(i) = .True.
    n_det_ref_tmp +=1
    index_ref_tmp(n_det_ref_tmp) = i
    coef_ref(n_det_ref_tmp) = psi_coef(i,1)
    exit
   endif
  enddo
 enddo
 if( ndet_generators_input .ne. n_det_ref_tmp)then
  print*,'Problem !!!! '
  print*,' ndet_generators .ne. n_det_ref_tmp !!!'
  print*,'ndet_generators,n_det_ref_tmp'
  print*,ndet_generators_input,n_det_ref_tmp 
  stop
 endif
 
 call i_h_j(psi_det_generators_input(1,1,1),psi_det_generators_input(1,1,1),N_int,href)
 integer :: i_pert, i_pert_count
 i_pert_count = 0
 do i = 1, N_det
  if(is_a_ref_det(i))cycle
  call i_h_j(psi_det(1,1,i),psi_det(1,1,i),N_int,hka)
  double precision :: f,href
  f = 1.d0/(href - hka)
  H_array = 0.d0
  accu = 0.d0
  do j=1,ndet_generators_input
   call i_h_j(psi_det(1,1,i),psi_det_generators_input(1,1,j),N_int,hka)
   H_array(j) = hka
   accu += coef_ref(j) * hka
  enddo
  lambda_i = psi_coef(i,1)/accu
  i_pert = 1
  if(accu * f / psi_coef(i,1) .gt. 0.5d0 .and. accu * f/psi_coef(i,1).gt.0.d0)then
   i_pert = 0
  endif
  do j = 1, ndet_generators_input
   if(dabs(H_array(j)*lambda_i).gt.0.1d0)then
    i_pert = 1
    exit
   endif
  enddo
  if(i_pert==1)then
   lambda_i = f
   i_pert_count +=1
  endif
  do k=1,ndet_generators_input
    double precision :: contrib
    contrib = H_array(k) * H_array(k) * lambda_i
    dressing_matrix(k, k) += contrib
    do j=k+1,ndet_generators_input
      contrib = H_array(k) * H_array(j) * lambda_i
      dressing_matrix(k, j) += contrib
      dressing_matrix(j, k) += contrib
    enddo 
  enddo
 enddo
 href = dressing_matrix(1,1)
 print*,'Diagonal part of the dressing'
 do i = 1, ndet_generators_input
  print*,'delta e = ',dressing_matrix(i,i) - href
 enddo
!print*,'i_pert_count = ',i_pert_count
end



