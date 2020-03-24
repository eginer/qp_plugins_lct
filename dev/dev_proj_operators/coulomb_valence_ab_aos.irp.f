BEGIN_PROVIDER [double precision, rho_alpha_hf_ao, (ao_num,ao_num)]
 implicit none 
 integer :: i,j,k,ii
 rho_alpha_hf_ao = 0.d0
 do i = 1, n_occ_val_orb_for_hf(1)
  ii = list_valence_orb_for_hf(i,1)
  do j = 1, ao_num
   do k = 1, ao_num
    rho_alpha_hf_ao(j,k) += mo_coef(j,ii) * mo_coef(k,ii)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, rho_alpha_hf_ao_in_r, (n_points_final_grid,ao_num) ]
 implicit none
 integer :: ipoint, j,k
 rho_alpha_hf_ao_in_r = 0.d0
 do j = 1, ao_num
  do ipoint = 1, n_points_final_grid
   do k = 1, ao_num
    rho_alpha_hf_ao_in_r(ipoint,j) += rho_alpha_hf_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, rho_alpha_hf_ao_in_r_bis, (ao_num,n_points_final_grid) ]
 implicit none
 integer :: ipoint, j,k
 rho_alpha_hf_ao_in_r_bis = 0.d0
 do ipoint = 1, n_points_final_grid
  do j = 1, ao_num
   do k = 1, ao_num
    rho_alpha_hf_ao_in_r_bis(j,ipoint) += rho_alpha_hf_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, rho_beta_hf_ao, (ao_num,ao_num)]
 implicit none 
 integer :: i,j,k,ii
 rho_beta_hf_ao = 0.d0
 do i = 1, n_occ_val_orb_for_hf(2)
  ii = list_valence_orb_for_hf(i,2)
  do j = 1, ao_num
   do k = 1, ao_num
    rho_beta_hf_ao(j,k) += mo_coef(j,ii) * mo_coef(k,ii)
   enddo
  enddo
 enddo

END_PROVIDER 

BEGIN_PROVIDER [double precision, rho_beta_hf_ao_in_r, (n_points_final_grid,ao_num) ]
 implicit none
 integer :: ipoint, j,k
 rho_beta_hf_ao_in_r = 0.d0
 do j = 1, ao_num
  do ipoint = 1, n_points_final_grid
   do k = 1, ao_num
    rho_beta_hf_ao_in_r(ipoint,j) += rho_beta_hf_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, rho_beta_hf_ao_in_r_bis, (ao_num,n_points_final_grid) ]
 implicit none
 integer :: ipoint, j,k
 rho_beta_hf_ao_in_r_bis = 0.d0
 do ipoint = 1, n_points_final_grid
  do j = 1, ao_num
   do k = 1, ao_num
    rho_beta_hf_ao_in_r_bis(j,ipoint) += rho_beta_hf_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, full_dens_ao, (ao_num,ao_num)]
 implicit none
 integer :: i,j,k 
 full_dens_ao = 0.d0
 do i = 1, mo_num
  do j = 1, ao_num
   do k = 1, ao_num
    full_dens_ao(k,j) += mo_coef(k,i) * mo_coef(j,i) 
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, full_dens_ao_in_r, (n_points_final_grid,ao_num) ]
 implicit none
 integer :: ipoint, j,k
 full_dens_ao_in_r = 0.d0
 do j = 1, ao_num
  do ipoint = 1, n_points_final_grid
   do k = 1, ao_num
    full_dens_ao_in_r(ipoint,j) += full_dens_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, full_dens_ao_in_r_bis, (ao_num,n_points_final_grid) ]
 implicit none
 integer :: ipoint, j,k
 full_dens_ao_in_r_bis = 0.d0
 do ipoint = 1, n_points_final_grid
  do j = 1, ao_num
   do k = 1, ao_num
    full_dens_ao_in_r_bis(j,ipoint) += full_dens_ao(k,j) * aos_in_r_array(k,ipoint)
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, f_hf_ab_ao, (n_points_final_grid)]
 implicit none
 integer :: alph,bet,gam,delt,i,ipoint,jpoint
 double precision :: thresh, integral 
 integer :: sze,sze_max,non_zero_int
 integer, allocatable :: out_val_index(:,:)
 double precision, allocatable :: out_val(:),alpha_dens(:)
 sze = ao_num
 sze_max = ao_num * ao_num
 thresh = thresh_int_mu_of_r

 f_hf_ab_ao = 0.d0
 provide rho_alpha_hf_ao_in_r full_dens_ao_in_r rho_beta_hf_ao_in_r ao_two_e_integrals_in_map
 print*,'all small stuff provided !'

 integer, allocatable  :: list_points(:)
 integer :: n_good_points
 allocate(list_points(n_points_final_grid))
 !  First alpha pair of AO
 do alph = 1, ao_num !
  do gam = 1, ao_num !
   if (ao_overlap_abs(alph,gam) < thresh) cycle

   allocate(out_val(sze_max),out_val_index(2,sze_max),alpha_dens(n_points_final_grid))
   call get_ao_two_e_integrals_non_zero_jl(alph,gam,thresh,sze_max,sze,out_val,out_val_index,non_zero_int) 
   n_good_points = 0
   do jpoint = 1, n_points_final_grid
    alpha_dens(jpoint) = rho_alpha_hf_ao_in_r(jpoint,gam) * full_dens_ao_in_r(jpoint,alph)
    if(dabs(alpha_dens(jpoint)).ge.thresh)then
     n_good_points += 1
     list_points(n_good_points) = jpoint 
    endif
   enddo
   print*,'n_good_points = ',n_good_points
   print*,'non_zero_int  = ',non_zero_int
   do i = 1, non_zero_int
    bet  = out_val_index(1,i)
    delt = out_val_index(2,i) 
    integral = out_val(i)

!!$OMP PARALLEL        &
!!$OMP DEFAULT (NONE)  &
!!$OMP SHARED (thresh,sze_max,sze,f_hf_ab_ao,n_points_final_grid,rho_alpha_hf_ao_in_r,full_dens_ao_in_r,rho_beta_hf_ao_in_r,ao_num)  & 
!!$OMP SHARED (alph,gam,bet,delt,integral,alpha_dens,out_val,out_val_index,non_zero_int) &
!!$OMP PRIVATE (ipoint)   
    ! OMP DO SCHEDULE (STATIC)
    do ipoint = 1, n_good_points
     f_hf_ab_ao(list_points(ipoint)) += rho_beta_hf_ao_in_r(list_points(ipoint),bet) * full_dens_ao_in_r(list_points(ipoint),delt) * integral * alpha_dens(list_points(ipoint)) 
    enddo
!   ! OMP END DO NO WAIT 
!!$OMP END PARALLEL
   enddo
   deallocate(out_val,out_val_index,alpha_dens)
  enddo
 enddo


 print*,'f_hf_ab_ao provided '

END_PROVIDER 


BEGIN_PROVIDER [double precision, f_hf_ab_ao_bis, (n_points_final_grid)]
 implicit none
 integer :: alph,bet,gam,delt,i,ipoint,jpoint
 double precision :: thresh, integral 
 integer :: sze,sze_max,non_zero_int,m
 integer, allocatable :: out_val_index(:)
 double precision, allocatable :: out_val(:),alpha_dens(:)
 sze = ao_num
 sze_max = ao_num
 thresh = thresh_int_mu_of_r

 f_hf_ab_ao_bis = 0.d0
 provide rho_alpha_hf_ao_in_r_bis full_dens_ao_in_r_bis rho_beta_hf_ao_in_r_bis ao_two_e_integrals_in_map
 print*,'all small stuff provided !'

 double precision :: tmp1,tmp2,tmp3,accu1,accu2,accu3
 allocate(out_val(sze_max),out_val_index(sze_max),alpha_dens(n_points_final_grid))
 !  First alpha pair of AO
 do ipoint = 1, n_points_final_grid
  do alph = 1, ao_num !
   tmp1 = dabs(rho_alpha_hf_ao_in_r_bis(alph,ipoint))
   accu1 = rho_alpha_hf_ao_in_r_bis(alph,ipoint)
   if( tmp1 .lt. thresh)cycle
   do gam = 1, ao_num !
    if (ao_overlap_abs(alph,gam) < thresh) cycle
    tmp2 = tmp1 * dabs(full_dens_ao_in_r_bis(gam,ipoint))
    accu2 = accu1 * full_dens_ao_in_r_bis(gam,ipoint)
    if (tmp2  .lt.thresh)cycle
    do bet = 1, ao_num
     tmp3 = tmp2 * dabs(rho_beta_hf_ao_in_r_bis(bet,ipoint))
     accu3 = accu2 * rho_beta_hf_ao_in_r_bis(bet,ipoint)
     if(tmp3 .lt. thresh)cycle
     call get_ao_two_e_integrals_non_zero(alph,bet,gam,sze,out_val,out_val_index,non_zero_int)
     do m = 1, non_zero_int
      delt = out_val_index(m) 
      integral = out_val(m) 
      f_hf_ab_ao_bis(ipoint) += accu3 * full_dens_ao_in_r_bis(delt,ipoint) 
     enddo
    enddo
   enddo
  enddo
 enddo


 print*,'f_hf_ab_ao_bis provided '

END_PROVIDER 



