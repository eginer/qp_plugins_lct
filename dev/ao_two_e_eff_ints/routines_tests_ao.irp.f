
subroutine test_gauss_ints_aos
 implicit none 
 double precision :: weight,r(3),weightj,r2(3),int_ao
 double precision :: ao_two_e_integral_schwartz_accel_gauss
 double precision :: int_r2,r_12,int_gauss_num,alpha_r12,coef
 integer :: ipoint,i,j,n_pt_in,jpoint
 double precision :: aos_array_r1(ao_num),aos_array_r2(ao_num)
 double precision :: accu_relat,accu_abs,err_relat,err_abs
 include 'utils/constants.include.F'
 integer :: iao,jao,kao,lao
 do iao = 1, ao_num ! r1
  do jao = 1, ao_num ! r2
   do kao = 1, ao_num ! r1
    do lao = 1, ao_num ! r2
 

     print*,'<ij|kl> = ',iao,jao,kao,lao
     int_ao = ao_two_e_integral_schwartz_accel_gauss(iao,kao,jao,lao)
     print*,'int_ao        = ',int_ao
     int_gauss_num = 0.d0
     do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      call give_all_aos_at_r(r,aos_array_r1)
      weight = final_weight_at_r_vector(ipoint)
      int_r2 = 0.d0
      do jpoint = 1, n_points_final_grid
       r2(1) = final_grid_points(1,jpoint)
       r2(2) = final_grid_points(2,jpoint)
       r2(3) = final_grid_points(3,jpoint)
       call give_all_aos_at_r(r2,aos_array_r2)
       weightj = final_weight_at_r_vector(jpoint)
       r_12 = (r(1) - r2(1))**2 + (r(2) - r2(2))**2 + (r(3) - r2(3))**2 
       do i = 1,n_gauss_eff_pot
        alpha_r12 = expo_gauss_eff_pot(i)
        if(alpha_r12 * r_12.gt.20.d0)cycle
        coef      = coef_gauss_eff_pot(i)
        int_r2   += aos_array_r2(jao) * aos_array_r2(lao) * dexp(-alpha_r12*r_12) * coef * weightj
       enddo
      enddo
      int_gauss_num += weight * int_r2 * aos_array_r1(iao) * aos_array_r1(kao)
     enddo
     err_abs = dabs(int_gauss_num - int_ao)
     if(int_gauss_num.gt.1.d-10)then
      err_relat = err_abs/dabs(int_gauss_num)
     else
      err_relat = 0.d0
     endif
     print*,'int_gauss_num = ',int_gauss_num
     print*,'abs error     = ',err_abs
     print*,'err_relat     = ',err_relat
     accu_abs += err_abs
     accu_relat += accu_relat
    enddo
   enddo
  enddo
 enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num**4)
 print*,'accu_relat = ',accu_relat/dble(ao_num**4)

end


subroutine test_extra_basis
 implicit none
 integer :: ipoint,i_ao,m
 double precision :: r(3),accu(3,ao_num),weight,aos_array(ao_num)
 double precision :: xyz_phi(3),aos_grad_array(3,ao_num),grad_xyz_phi(3),xyz_grad_phi(3)
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)
  weight = final_weight_at_r_vector(ipoint)
  do i_ao = 1, ao_num
   call xyz_phi_ao(r,i_ao,xyz_phi)
   call xyz_grad_phi_ao(r,i_ao,xyz_grad_phi)
   do m = 1, 3
!    if(dabs(aos_array(i_ao) * r(m) - xyz_phi(m)).gt.1.d-10)then
!    if(dabs(aos_grad_array(m,i_ao) - grad_xyz_phi(m)).gt.1.d-10)then
    if(dabs(r(m) * aos_grad_array(m,i_ao) - xyz_grad_phi(m)).gt.1.d-10)then
!     print*,i_ao,m,r(m)
     print*,i_ao,m,ao_power(i_ao,m)
     print*,r
     print*,aos_grad_array(m,i_ao)*r(m),xyz_grad_phi(m)
    endif
    accu(m,i_ao) += dabs(aos_grad_array(m,i_ao) * r(m) - xyz_grad_phi(m)) * weight
   enddo
  enddo
 enddo
 print*,''
 print*,'errors '
 print*,''
 do i_ao = 1, ao_num
  write(*,'(100(F16.10,X))')accu(:,i_ao)
 enddo
end

subroutine test_extra
 implicit none
 double precision :: r(3)
 integer :: prim_num_i,prim_num_j
 include 'utils/constants.include.F'
 double precision, allocatable :: P_new(:,:,:,:) ! new polynom for each couple of prim
 double precision :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 integer :: i,j,ipoint,num_i,num_j
 double precision :: ao_prod_in_r,weight,accu,num,ref
 double precision :: aos_grad_array(3,ao_num),aos_array(ao_num)
 allocate(P_new(0:max_dim,3,ao_prim_num_max,ao_prim_num_max)) ! new polynom for each couple of prim

 do i = 1, ao_num
  do j = 1, ao_num
   print*,'i,j',i,j
   call give_poly_ij(i,j,P_new,P_center,p_exp,fact_p,iorder_p,coef_prod)
   accu = 0.d0
   do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    call give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)
    num = ao_prod_in_r(r,ao_prim_num(i),ao_prim_num(j),P_new,P_center,p_exp,fact_p,iorder_p,coef_prod)
    ref = aos_array(i) * aos_array(j)
    if(dabs(num - ref).gt.1.d-9)then
     print*,r
     num_i = ao_nucl(i)
     num_j = ao_nucl(j)
     print*,'num_i,num_j',num_i,num_j
     print*,'ao_power i ', ao_power(i,:)
     print*,'ao_power j ', ao_power(j,:)
     print*,num,ref
     print*,'STOOOOOP '
     stop
    endif
    accu += dabs(num - ref) * weight
   enddo
   print*,'accu = ',accu
   if(dabs(accu).gt.1.d-9)then
    num_i = ao_nucl(i)
    num_j = ao_nucl(j)
    print*,'num_i,num_j',num_i,num_j
    print*,'ao_power i ', ao_power(i,:)
    print*,'ao_power j ', ao_power(j,:)
    print*,'STOOOOOP '
    stop
   endif
  enddo
 enddo



end
