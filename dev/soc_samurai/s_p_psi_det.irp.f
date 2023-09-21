  use bitmasks ! you need to include the bitmasks_module.f90 features

 BEGIN_PROVIDER [ integer, n_det_per_conf_s_p, (N_configuration)]
&BEGIN_PROVIDER [ integer, n_det_max_per_conf_s_p]
&BEGIN_PROVIDER [ integer, size_n_det_max_per_conf_sp ]
 implicit none
 BEGIN_DOC
! number of determinants for each configuration when S^+ is applied on psi_det
 END_DOC
 integer :: i, ndet_conf_sp,n_alpha_p
 integer(bit_kind) :: det_tmp(N_int,2)
 do i = 1, N_configuration
  det_tmp(:,:) = psi_configuration(:,:,i)
  n_alpha_p = elec_alpha_num + 1 
  call configuration_to_dets_size(psi_configuration(1,1,i),ndet_conf_sp,n_alpha_p,N_int)
  n_det_per_conf_s_p(i) = ndet_conf_sp
 enddo
 n_det_max_per_conf_s_p = maxval(n_det_per_conf_s_p)
 size_n_det_max_per_conf_sp = max(1,n_det_max_per_conf_s_p)

END_PROVIDER 

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_s_p, (N_int,2,size_n_det_max_per_conf_sp,N_configuration)]
&BEGIN_PROVIDER [ double precision,  psi_coef_s_p, (size_n_det_max_per_conf_sp,N_configuration,n_states)]
&BEGIN_PROVIDER [ double precision,  psi_norm_s_p, (n_states) ]
 implicit none
 BEGIN_DOC
 ! psi_det_s_p(:,:,j,i) = jth determinant belonging to the ith configuration when S^+ is applied on the whole list of 
 ! Slater determinants in psi_det
 !
 ! psi_coef_s_p(j,i,k) = coefficient of the jth determinant belonging to the ith configuration when S^+ is applied on the kth state 
 !
 ! psi_norm_s_p(k) = norm of the k-th state when S^+ is applied on the k-th state (should be 1 or 0)
 END_DOC
 double precision :: coef_sp(N_states,size_n_det_max_per_conf_sp),phase(elec_beta_num)
 integer(bit_kind) :: dets_sp(N_int,2,size_n_det_max_per_conf_sp),det_out(N_int,2,elec_beta_num),det_tmp(N_int,2)
 integer :: i,j,k,l,m,n_alpha_p,ndet_out,istate,degree,ndet_conf_sp
 double precision :: S,ms,factor_s_p,coef_tmp
 ms = 0.5d0 * dble(elec_alpha_num - elec_beta_num)
 do i = 1, N_configuration ! first loop over configurations (because S^+ does not couple different configurations)
  det_tmp(:,:) = psi_configuration(:,:,i)
  coef_sp = 0.d0
  n_alpha_p = elec_alpha_num + 1 ! Applying S^+ I increase by one the number of alpha electron
  ndet_conf_sp = n_det_per_conf_s_p(i) ! number of determinants in configuration i
  if(ndet_conf_sp  == 0)then
   psi_det_s_p(:,:,1,i) = 0_bit_kind
   psi_coef_s_p(1,i,:) = 0.d0
  else
   ! this routine returns all possible determinants belonging to the configuration i with n_alpha_p alpha electrons 
   call configuration_to_dets(psi_configuration(1,1,i),dets_sp,ndet_conf_sp,n_alpha_p,N_int)
   ! the set of determinants psi_det_s_p is the basis for S^+ on the whole set of determinants belonging to conf i
   do j = 1, ndet_conf_sp ! number of determinant corresponding to the configuration i with n_alpha +=1
    psi_det_s_p(:,:,j,i) = dets_sp(:,:,j)
   enddo
   ! loop over all the determinants in psi_det belonging to configuration i 
   do j =  psi_configuration_to_psi_det(1,i), psi_configuration_to_psi_det(2,i)
    k = psi_configuration_to_psi_det_data(j) ! index of the determinant "j" in the configuration i 
    call s_plus_det(psi_det(1,1,k),det_out,phase,ndet_out) ! S^+ |j> = \sum_l=1,ndet_out phase_l |det_out_l>
    ! Then you increment the coefficients with the phase in ndet_conf_sp
    do m = 1, ndet_conf_sp ! loop over the basis of the conf i with n_alpha+=1
     do l = 1, ndet_out ! loop over the S^+ |j> 
      call get_excitation_degree(dets_sp(1,1,m),det_out(1,1,l),degree,N_int)
      if(degree == 0)then ! if |det_out_l> = | det(i)> then 
       do istate = 1, N_states
        coef_sp(istate,m) += psi_coef(k,istate) * phase(m) ! increment the initial coefficient multiplied by the phase
       enddo
      endif
     enddo
    enddo
   enddo
   ! I divide by the prefactor sqrt(S(S+1) - S_z(S_z+1)) in order to have a normalized state 
   do j = 1, n_det_per_conf_s_p(i) ! 
    do istate = 1, N_states
     S = s_values(istate)
     coef_tmp = factor_s_p(S,ms)
     if(dabs(coef_tmp).gt.1.d-12)then
      psi_coef_s_p(j,i,istate) = coef_sp(istate,j) / coef_tmp
     else
      psi_coef_s_p(j,i,istate) = coef_sp(istate,j) 
     endif
    enddo
   enddo
  endif
 enddo
 psi_norm_s_p = 0.d0
 do istate = 1, N_states
  do i = 1, N_configuration
   do j = 1, n_det_per_conf_s_p(i)
    psi_norm_s_p(istate) += psi_coef_s_p(j,i,istate)**2.d0
   enddo
  enddo
 enddo

 END_PROVIDER 
