
  use bitmasks ! you need to include the bitmasks_module.f90 features

 BEGIN_PROVIDER [ integer, n_det_per_conf_s_m, (N_configuration)]
&BEGIN_PROVIDER [ integer, n_det_max_per_conf_s_m]
&BEGIN_PROVIDER [ integer, size_n_det_max_per_conf_sm ]
 implicit none
 BEGIN_DOC
! number of determinants for each configuration when S^- is applied on psi_det
 END_DOC
 integer :: i, ndet_conf_sm,n_alpha_m
 integer(bit_kind) :: det_tmp(N_int,2)
 do i = 1, N_configuration
  det_tmp(:,:) = psi_configuration(:,:,i)
  n_alpha_m = elec_alpha_num - 1 
  n_alpha_m = max(1,n_alpha_m)
  call configuration_to_dets_size(psi_configuration(1,1,i),ndet_conf_sm,n_alpha_m,N_int)
  n_det_per_conf_s_m(i) = ndet_conf_sm
  print*,'n_det_per_conf_s_m = ',n_det_per_conf_s_m(i)
 enddo
 n_det_max_per_conf_s_m = maxval(n_det_per_conf_s_m)
 size_n_det_max_per_conf_sm = max(1,n_det_max_per_conf_s_m)

END_PROVIDER 

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_s_m, (N_int,2,size_n_det_max_per_conf_sm,N_configuration)]
&BEGIN_PROVIDER [ double precision,  psi_coef_s_m, (size_n_det_max_per_conf_sm,N_configuration,n_states)]
&BEGIN_PROVIDER [ double precision,  psi_norm_s_m, (n_states) ]
 implicit none
 double precision :: coef_sm(N_states,size_n_det_max_per_conf_sm),phase(elec_alpha_num)
 integer(bit_kind) :: dets_sm(N_int,2,size_n_det_max_per_conf_sm),det_out(N_int,2,elec_alpha_num),det_tmp(N_int,2)
 integer :: i,j,k,l,m,n_alpha_m,ndet_out,istate,degree,ndet_conf_sm
 double precision :: S,ms,factor_s_m,coef_tmp
 ms = 0.5d0 * dble(elec_alpha_num - elec_beta_num)
 do i = 1, N_configuration ! first loop over configurations (because S^- does not couple different configurations)
  det_tmp(:,:) = psi_configuration(:,:,i)
  coef_sm = 0.d0
  n_alpha_m = elec_alpha_num - 1 ! Applying S^- I decrease by one the number of alpha electron
  if(elec_alpha_num == 1)then ! just have to spin-flip the determinant 
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
  else
   ndet_conf_sm = n_det_per_conf_s_m(i)
   if(ndet_conf_sm == 0)then
    psi_det_s_m(:,:,1,i) = 0_bit_kind
    psi_coef_s_m(1,i,:) = 0.d0
   else
    if(elec_alpha_num .gt. 1 .and. elec_beta_num .gt. 0)then
     call configuration_to_dets(psi_configuration(1,1,i),dets_sm,ndet_conf_sm,n_alpha_m,N_int)
     do j = 1, ndet_conf_sm
      psi_det_s_m(:,:,j,i) = dets_sm(:,:,j)
     enddo
    else if(elec_alpha_num == 1. and. elec_alpha_num == 0)then ! special case for elec_alpha_num == 1
     do j =  psi_configuration_to_psi_det(1,i), psi_configuration_to_psi_det(2,i)
      k = psi_configuration_to_psi_det_data(j)
      call s_minus_det(psi_det(1,1,k),det_out,phase,ndet_out)
      psi_det_s_m(:,:,j,i) = det_out(:,:)
     enddo 
    endif
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
   endif
  endif
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
