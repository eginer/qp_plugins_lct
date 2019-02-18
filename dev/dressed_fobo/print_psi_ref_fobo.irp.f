program pouet
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end
 
subroutine routine
 implicit none
 integer :: i,j
 integer :: number_of_holes,nh
 integer :: number_of_particles,np
 do i =1, N_det
  nh = number_of_holes(psi_det(1,1,i))
  np = number_of_particles(psi_det(1,1,i))
  double precision :: psi_coef_av
  psi_coef_av = 0.d0
  if(nh == 2 .and. np == 0 )then
   print*,'2h  ',psi_coef(i,:)
  endif
  do j = 1, N_states
   psi_coef_av += dabs(psi_coef(i,j))
  enddo
  psi_coef_av = psi_coef_av / dble(N_states)
  if(psi_coef_av.lt.thresh_coef_lmct)cycle
  if(     nh == 0 .and. np == 0 )then
   print*,'CAS ', psi_coef(i,:) 
  else if(nh == 1 .and. np == 0 )then
   print*,'1h  ',psi_coef(i,:)
  endif
 enddo


end

subroutine routine_2
 implicit none
 integer :: i
 print*, 'N_det_ref_fobo     = ',N_det_ref_fobo
 print*, 'N_det_non_ref_fobo = ',N_det_non_ref_fobo
 print*, 'N_det_lmct         = ',N_det_lmct
 print*, 'N_det_mlct         = ',N_det_mlct
 print*, '-----'
 print*, '-----'
 print*, '     '
 print*, 'psi_det fobo'
 do i = 1, N_det_ref_fobo
 print*, '-----'
  print*, 'i',i
  print*, 'coef = ',psi_ref_fobo_coef(i,:)
  call debug_det(psi_ref_fobo(1,1,i),N_int)
 print*, '-----'
 enddo

end
