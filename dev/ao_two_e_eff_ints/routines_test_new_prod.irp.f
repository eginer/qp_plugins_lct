
subroutine test_all_poly_for_r12_deriv(i,j,r,value_ij,value_j_x_i, value_j_dxyz_i, value_j_xyz_dxyz_i)
 include 'utils/constants.include.F'                                                                                                                                  
 implicit none
 integer, intent(in) :: i,j
 double precision, intent(in) :: r(3)
 double precision, intent(out):: value_ij,value_j_x_i(3),value_j_dxyz_i(3),value_j_xyz_dxyz_i(3)

 double precision :: P_ij(0:max_dim,3,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer          :: iorder_ij(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision :: P_center_ij(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision :: p_exp_ij(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision :: fact_p_ij(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 double precision :: coef_prod_ij(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision :: P_j_xyz_i(0:max_dim,3,2,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer          :: iorder_j_xyz_i(3,2,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision :: P_j_dxyz_i(0:max_dim,3,2,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer          :: iorder_j_dxyz_i(3,2,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision :: P_j_xyz_dxyz_i(0:max_dim,3,4,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer          :: iorder_j_xyz_dxyz_i(3,4,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim

 call give_all_poly_for_r12_deriv(i,j,P_ij,iorder_ij,P_center_ij,p_exp_ij,fact_p_ij,coef_prod_ij,& 
            P_j_xyz_i,iorder_j_xyz_i, iorder_j_dxyz_i,P_j_dxyz_i, P_j_xyz_dxyz_i, iorder_j_xyz_dxyz_i)

 double precision :: ao_prod_in_r_bis
 double precision :: P_x(0:max_dim,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 double precision :: P_y(0:max_dim,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim                                                  
 double precision :: P_z(0:max_dim,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer          :: iorder_tmp(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 
 integer :: p,q,prim_num_i,prim_num_j
 prim_num_i = ao_prim_num(i)
 prim_num_j = ao_prim_num(j)
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_ij(0:max_dim,1,q,p)
   P_y(0:max_dim,q,p) = P_ij(0:max_dim,2,q,p)
   P_z(0:max_dim,q,p) = P_ij(0:max_dim,3,q,p)
   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! phi_j *  phi_i = (phi_j * phi_i)_x * (phi_j * phi_i)_y * (phi_j * phi_i)_z 
 value_ij = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! phi_j * x * phi_i = (phi_j * phi_i)_y * (phi_j * phi_i)_z * (phi_j * x * phi_i)_x
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_j_xyz_i(0:max_dim,1,1,q,p) ! the part with the new polynoms in x 
   P_y(0:max_dim,q,p) = P_ij(0:max_dim,2,q,p) ! the usual part with the polynom in y 
   P_z(0:max_dim,q,p) = P_ij(0:max_dim,3,q,p) ! the usual part with the polynom in z
   iorder_tmp(1,q,p)  = iorder_j_xyz_i(1,1,q,p)
   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_x_i(1) = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_j_xyz_i(0:max_dim,1,2,q,p) ! you replace by the second contribution in x 
   iorder_tmp(1,q,p)  = iorder_j_xyz_i(1,2,q,p)
!   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
!   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_x_i(1) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)

 ! phi_j * y * phi_i = (phi_j * phi_i)_x * (phi_j * phi_i)_z * (phi_j * y * phi_i)_y
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_ij(0:max_dim,1,q,p) ! the usual part with the polynom in x 
   P_y(0:max_dim,q,p) = P_j_xyz_i(0:max_dim,2,1,q,p) ! the part with the new polynoms in y 
   P_z(0:max_dim,q,p) = P_ij(0:max_dim,3,q,p) ! the usual part with the polynom in z
   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
   iorder_tmp(2,q,p)  = iorder_j_xyz_i(2,1,q,p)
   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_x_i(2) = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_y(0:max_dim,q,p) = P_j_xyz_i(0:max_dim,2,2,q,p) ! you replace by the second contribution in y 
!   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
   iorder_tmp(2,q,p)  = iorder_j_xyz_i(2,2,q,p)
!   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_x_i(2) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)

 ! phi_j * z * phi_i = (phi_j * phi_i)_x * (phi_j * phi_i)_y * (phi_j * z * phi_i)_z
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_ij(0:max_dim,1,q,p) ! the usual part with the polynom in x 
   P_y(0:max_dim,q,p) = P_ij(0:max_dim,2,q,p) ! the usual part with the polynom in y
   P_z(0:max_dim,q,p) = P_j_xyz_i(0:max_dim,3,1,q,p) ! the part with the new polynoms in z 
   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
   iorder_tmp(3,q,p)  = iorder_j_xyz_i(3,1,q,p)
  enddo
 enddo
 value_j_x_i(3) = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_z(0:max_dim,q,p) = P_j_xyz_i(0:max_dim,3,2,q,p) ! you replace by the second contribution in z 
!   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
!   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
   iorder_tmp(3,q,p)  = iorder_j_xyz_i(3,2,q,p)
  enddo
 enddo
 value_j_x_i(3) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! phi_j * d/dx * phi_i = (phi_j * phi_i)_y * (phi_j * phi_i)_z * (phi_j * d/dx * phi_i)_x
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_j_dxyz_i(0:max_dim,1,1,q,p) ! the part with the new polynoms in x 
   P_y(0:max_dim,q,p) = P_ij(0:max_dim,2,q,p) ! the usual part with the polynom in y 
   P_z(0:max_dim,q,p) = P_ij(0:max_dim,3,q,p) ! the usual part with the polynom in z
   iorder_tmp(1,q,p)  = iorder_j_dxyz_i(1,1,q,p)
   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_dxyz_i(1) = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_j_dxyz_i(0:max_dim,1,2,q,p) ! you replace by the second contribution in x 
   iorder_tmp(1,q,p)  = iorder_j_dxyz_i(1,2,q,p)
!   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
!   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_dxyz_i(1) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)

 ! phi_j * d/dy * phi_i = (phi_j * phi_i)_x * (phi_j * phi_i)_z * (phi_j * d/dy * phi_i)_y
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_ij(0:max_dim,1,q,p) ! the usual part with the polynom in x 
   P_y(0:max_dim,q,p) = P_j_dxyz_i(0:max_dim,2,1,q,p) ! the part with the new polynoms in y 
   P_z(0:max_dim,q,p) = P_ij(0:max_dim,3,q,p) ! the usual part with the polynom in z
   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
   iorder_tmp(2,q,p)  = iorder_j_dxyz_i(2,1,q,p)
   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_dxyz_i(2) = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_y(0:max_dim,q,p) = P_j_dxyz_i(0:max_dim,2,2,q,p) ! you replace by the second contribution in y 
!   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
   iorder_tmp(2,q,p)  = iorder_j_dxyz_i(2,2,q,p)
!   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
  enddo
 enddo
 value_j_dxyz_i(2) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)

 ! phi_j * d/dz * phi_i = (phi_j * phi_i)_x * (phi_j * phi_i)_y * (phi_j * d/dz * phi_i)_z
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_x(0:max_dim,q,p) = P_ij(0:max_dim,1,q,p) ! the usual part with the polynom in x 
   P_y(0:max_dim,q,p) = P_ij(0:max_dim,2,q,p) ! the usual part with the polynom in y
   P_z(0:max_dim,q,p) = P_j_dxyz_i(0:max_dim,3,1,q,p) ! the part with the new polynoms in z 
   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
   iorder_tmp(3,q,p)  = iorder_j_dxyz_i(3,1,q,p)
  enddo
 enddo
 value_j_dxyz_i(3) = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
 do p = 1, prim_num_j 
  do q = 1, prim_num_i 
   P_z(0:max_dim,q,p) = P_j_dxyz_i(0:max_dim,3,2,q,p) ! you replace by the second contribution in z 
!   iorder_tmp(1,q,p)  = iorder_ij(1,q,p)
!   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
   iorder_tmp(3,q,p)  = iorder_j_dxyz_i(3,2,q,p)
  enddo
 enddo
 value_j_dxyz_i(3) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! phi_j * x * d/dx * phi_i = (phi_j * phi_i)_y * (phi_j * phi_i)_z * (phi_j * x * d/dx * phi_i)_x
! do p = 1, prim_num_j 
!  do q = 1, prim_num_i 
!   P_x(0:max_dim,q,p) = P_j_xyz_dxyz_i(0:max_dim,1,1,q,p) ! the part with the new polynoms in x 
!   P_y(0:max_dim,q,p) = P_ij(0:max_dim,2,q,p) ! the usual part with the polynom in y 
!   P_z(0:max_dim,q,p) = P_ij(0:max_dim,3,q,p) ! the usual part with the polynom in z
!   iorder_tmp(1,q,p)  = iorder_j_xyz_dxyz_i(1,1,q,p)
!   iorder_tmp(2,q,p)  = iorder_ij(2,q,p)
!   iorder_tmp(3,q,p)  = iorder_ij(3,q,p)
!  enddo
! enddo
! value_j_xyz_dxyz_i(1) = ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
! do p = 1, prim_num_j 
!  do q = 1, prim_num_i 
!   P_x(0:max_dim,q,p) = P_j_xyz_dxyz_i(0:max_dim,1,3,q,p) ! you replace by the second contribution in x 
!   iorder_tmp(1,q,p)  = iorder_j_xyz_dxyz_i(1,3,q,p)
!  enddo
! enddo
! value_j_xyz_dxyz_i(1) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
! do p = 1, prim_num_j 
!  do q = 1, prim_num_i 
!   P_x(0:max_dim,q,p) = P_j_xyz_dxyz_i(0:max_dim,1,4,q,p) ! you replace by the second contribution in x 
!   iorder_tmp(1,q,p)  = iorder_j_xyz_dxyz_i(1,4,q,p)
!  enddo
! enddo
! value_j_xyz_dxyz_i(1) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)

 integer :: kk
 value_j_xyz_dxyz_i = 0.d0
 do kk = 1, 4 ! loop on the number of different polynoms 
  do p = 1, prim_num_j 
   do q = 1, prim_num_i 
    P_x(0:max_dim,q,p) = P_j_xyz_dxyz_i(0:max_dim,1,kk,q,p) ! you replace by the second contribution in x 
    P_y(0:max_dim,q,p) = P_ij(0:max_dim,2,q,p) ! the usual part with the polynom in y 
    P_z(0:max_dim,q,p) = P_ij(0:max_dim,3,q,p) ! the usual part with the polynom in z
    iorder_tmp(1,q,p)  = iorder_j_xyz_dxyz_i(1,kk,q,p) ! you replace the order of polynom by that in x 
    iorder_tmp(2,q,p)  = iorder_ij(2,q,p) ! usual polynom for prod i * j in y
    iorder_tmp(3,q,p)  = iorder_ij(3,q,p) ! usual polynom for prod i * j in z
   enddo
  enddo
  value_j_xyz_dxyz_i(1) += ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center_ij,p_exp_ij,fact_p_ij,iorder_tmp,coef_prod_ij)
 enddo


end


subroutine test_all_prod_in_r
 implicit none
 double precision :: r(3)
 integer :: prim_num_i,prim_num_j
 include 'utils/constants.include.F'
 double precision :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 integer :: i,j,ipoint,num_i,num_j,k
 double precision :: ao_prod_in_r,weight,accu(3),num,ref,accu_scal
 double precision :: aos_grad_array(3,ao_num),aos_array(ao_num),aos_dxyzi_j(3),phi_ao_plus_n,ao_i1
 double precision :: value_ij,value_j_x_i(3),value_j_dxyz_i(3),value_j_xyz_dxyz_i(3)

 do i = 1, ao_num
  do j = 1, ao_num
! do i = 1, ao_num
!  do j = 71, ao_num
   accu = 0.d0
   accu_scal = 0.d0
   do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    call give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)
    call test_all_poly_for_r12_deriv(i,j,r,value_ij,value_j_x_i,value_j_dxyz_i,value_j_xyz_dxyz_i)
!    ref = aos_array(i) * aos_array(j) 
!    num = value_ij
!    accu_scal += dabs(num - ref) * weight
    do k = 1,1
!     ref = aos_array(i) * aos_array(j) * r(k)
!     ref = aos_grad_array(k,i) * aos_array(j)
     ref = aos_grad_array(k,i) * aos_array(j) * r(k)
!     num = value_j_x_i(k)
!     num = value_j_dxyz_i(k)
     num = value_j_xyz_dxyz_i(k)
     if(dabs(num - ref).gt.1.d-9)then
      print*,'test_prod_xyz_dxyzi_j'
      print*,'STOOOOOP '
      print*,'STOOOOOP '
      print*,'i,j',i,j
      print*,r
      num_i = ao_nucl(i)
      num_j = ao_nucl(j)
      print*,'k = ',k
      print*,'num_i,num_j',num_i,num_j
      print*,'ao_power i ', ao_power(i,:)
      print*,'ao_power j ', ao_power(j,:)
      print*,num,ref
      print*,'STOOOOOP '
      print*,ipoint
      print*,'grad_i = ',aos_grad_array(k,i)
      print*,'ao j   = ',aos_array(j)
      stop
     endif
     accu(k) += dabs(num - ref) * weight
    enddo
   enddo
   print*,'i,j',i,j
!   print*,accu_scal
   print*,accu
  enddo
 enddo
end

