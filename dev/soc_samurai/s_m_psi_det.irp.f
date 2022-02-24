
  use bitmasks ! you need to include the bitmasks_module.f90 features

 BEGIN_PROVIDER [ integer, n_det_per_conf_s_m, (N_configuration)]
&BEGIN_PROVIDER [ integer, n_det_max_per_conf_s_m]
 implicit none
 BEGIN_DOC
! number of determinants for each configuration when S^- is applied on psi_det
 END_DOC
 integer :: i, ndet_conf_sm,n_alpha_p
 integer(bit_kind) :: det_tmp(N_int,2)
 do i = 1, N_configuration
  det_tmp(:,:) = psi_configuration(:,:,i)
  n_alpha_p = elec_alpha_num - 1 
  call configuration_to_dets_size(psi_configuration(1,1,i),ndet_conf_sm,n_alpha_p,N_int)
  n_det_per_conf_s_m(i) = ndet_conf_sm
 enddo
 n_det_max_per_conf_s_m = maxval(n_det_per_conf_s_m)

END_PROVIDER 

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_s_m, (N_int,2,n_det_max_per_conf_s_m,N_configuration)]
&BEGIN_PROVIDER [ double precision,  psi_coef_s_m, (n_det_max_per_conf_s_m,N_configuration,n_states)]
&BEGIN_PROVIDER [ double precision,  psi_norm_s_m, (n_states) ]
 implicit none
 double precision :: coef_sm(N_states,n_det_max_per_conf_s_m),phase(elec_beta_num)
 integer(bit_kind) :: dets_sm(N_int,2,n_det_max_per_conf_s_m),det_out(N_int,2,elec_beta_num),det_tmp(N_int,2)
 integer :: i,j,k,l,m,n_alpha_p,ndet_out,istate,degree,ndet_conf_sm
 double precision :: S,ms,factor_s_m,coef_tmp
 ms = 0.5d0 * dble(elec_alpha_num - elec_beta_num)
 do i = 1, N_configuration
  det_tmp(:,:) = psi_configuration(:,:,i)
  coef_sm = 0.d0
  n_alpha_p = elec_alpha_num - 1
  ndet_conf_sm = n_det_per_conf_s_m(i)
  call configuration_to_dets(psi_configuration(1,1,i),dets_sm,ndet_conf_sm,n_alpha_p,N_int)
  do j = 1, ndet_conf_sm
   psi_det_s_m(:,:,j,i) = dets_sm(:,:,j)
  enddo
  do j =  psi_configuration_to_psi_det(1,i), psi_configuration_to_psi_det(2,i)
   k = psi_configuration_to_psi_det_data(j)
   call s_minus_det(psi_det(1,1,k),det_out,phase,ndet_out)
   do m = 1, ndet_conf_sm
    do l = 1, ndet_out
     call get_excitation_degree(dets_sm(1,1,m),det_out(1,1,l),degree,N_int)
     if(degree == 0)then
      do istate = 1, N_states
       coef_sm(istate,m) += psi_coef(k,istate) * phase(m)
      enddo
     endif
    enddo
   enddo
  enddo
  do j = 1, n_det_per_conf_s_m(i)
   do istate = 1, N_states
    S = s_values(istate)
    coef_tmp = factor_s_m(S,ms)
    if(dabs(coef_tmp).gt.1.d-12)then
     psi_coef_s_m(j,i,istate) = coef_sm(istate,j) / coef_tmp
    else
     psi_coef_s_m(j,i,istate) = coef_sm(istate,j) 
    endif
   enddo
  enddo
 enddo
 psi_norm_s_m = 0.d0
 do istate = 1, N_states
  do i = 1, N_configuration
   do j = 1, n_det_per_conf_s_m(i)
    psi_norm_s_m(istate) += psi_coef_s_m(j,i,istate)**2.d0
   enddo
  enddo
 enddo

 END_PROVIDER 
