program soc_samurai
  implicit none
  BEGIN_DOC
  !TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  
  read_wf = .True.
  touch read_wf 
!  provide v_c_ij_grid 
!  provide D_mu_nu_find_a_better_name 
!  provide One_body_mu_nu_find_a_better_name 
end


subroutine test_sz
 implicit none
 integer :: i,j,k,ii,jj,kk
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer(bit_kind), allocatable :: det_i(:,:),det_j(:,:),conf_i(:,:), conf_j(:,:)
 double precision :: soc_ij
 allocate(det_i(N_int,2),conf_i(N_int,2))
 allocate(det_j(N_jnt,2),conf_j(N_jnt,2))
 ! Browse psi_det
 do i = 1, N_configuration ! loop over configurations :: same spatial part 
  conf_i(:,:) = psi_configuration(:,:,i) ! conf_i == spatial configuration "i"
  do j =  psi_configuration_to_psi_det(1,i), psi_configuration_to_psi_det(2,i) ! loop over determinants sharing the same spatial part 
   k = psi_configuration_to_psi_det_data(j) ! index of the determinant "j" in the configuration i 
   det_i(:,:) = psi_det(:,:,k)
   do ii = 1, N_configuration ! loop over configurations :: same spatial part 
    conf_j(:,:) = psi_configuration(:,:,ii) ! conf_i == spatial configuration "i"
    do jj =  psi_configuration_to_psi_det(1,ii), psi_configuration_to_psi_det(2,ii) ! loop over determinants sharing the same spatial part 
     kk = psi_configuration_to_psi_det_data(jj) ! index of the determinant "j" in the configuration i 
     det_j(:,:) = psi_det(:,:,kk)
     call i_soc_j(conf_i,conf_j,det_i, det_j,sz_i,sz_j,soc_ij)
    enddo
   enddo
  enddo
 enddo
end
