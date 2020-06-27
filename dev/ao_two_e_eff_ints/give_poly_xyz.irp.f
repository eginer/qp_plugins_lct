subroutine give_all_poly_for_r12_deriv(i,j,P_ij,iorder_ij,P_center_ij,p_exp_ij,fact_p_ij,coef_prod_ij,& 
            P_j_xyz_i,iorder_j_xyz_i, iorder_j_dxyz_i,P_j_dxyz_i, P_j_xyz_dxyz_i, iorder_j_xyz_dxyz_i)
 include 'utils/constants.include.F'                                                                                                                                  
 implicit none
 integer, intent(in) :: i,j
 double precision, intent(out) :: P_ij(0:max_dim,3,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer, intent(out)          :: iorder_ij(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: P_center_ij(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(out) :: p_exp_ij(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(out) :: fact_p_ij(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 double precision, intent(out) :: coef_prod_ij(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision, intent(out) :: P_j_xyz_i(0:max_dim,3,2,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer, intent(out)          :: iorder_j_xyz_i(3,2,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: P_j_dxyz_i(0:max_dim,3,2,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer, intent(out)          :: iorder_j_dxyz_i(3,2,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: P_j_xyz_dxyz_i(0:max_dim,3,4,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer, intent(out)          :: iorder_j_xyz_dxyz_i(3,4,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim

 integer :: p,q,k,m,n,num_ao,nx
 double precision :: center_ao(1:3)
 num_ao = ao_nucl(i)
 center_ao(1:3) = nucl_coord(num_ao,1:3)

 call give_poly_ij(i,j,P_ij,P_center_ij,p_exp_ij,fact_p_ij,iorder_ij,coef_prod_ij)

 ! phi_j x1 phi_i
 ! First contribution is A_x (x-A_x)^a_x
 do p = 1, ao_prim_num(j)
  do q = 1, ao_prim_num(i)
   do k = 1, 3
    iorder_j_xyz_i(k,1,q,p) = iorder_ij(k,q,p)
    do m = 0, iorder_ij(k,q,p)
     ! This corresponds to A_x (x - A_x)^a_x 
     P_j_xyz_i(m,k,1,q,p) = P_ij(m,k,q,p) * center_ao(k) 
    enddo
   enddo
  enddo
 enddo

 ! Second contribution is (x-A_x)^{a_x+1}
 integer :: n_new(3)
 n_new = 1
 double precision :: P_new(0:max_dim,3,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 ! You increase the order of the polynom of phi_i on x/y/z
 ! and you create all polynoms for the x/y/z
 call give_poly_i_plus_n_j_xyz(i,j,P_new,iorder_p,n_new)
 do p = 1, ao_prim_num(j)
  do q = 1, ao_prim_num(i)
   do k = 1, 3
    iorder_j_xyz_i(k,2,q,p) = iorder_p(k,q,p)
    do m = 0, iorder_p(k,q,p)
     ! This corresponds to (x-A_x)^{a_x+1}
     P_j_xyz_i(m,k,2,q,p) = P_new(m,k,q,p) 
    enddo
   enddo
  enddo
 enddo

 ! phi_j d/dx phi_i
 ! First contribution is a_x (x - A_x)^{a_x-1}
 n_new = -1
 call give_poly_i_plus_n_j_xyz(i,j,P_new,iorder_p,n_new)
 do p = 1, ao_prim_num(j)
  do q = 1, ao_prim_num(i)
   do k = 1, 3
    nx = ao_power(i,k)
    iorder_j_dxyz_i(k,1,q,p) = iorder_p(k,q,p)
    do m = 0, iorder_p(k,q,p)
     ! This corresponds to (x-A_x)^{a_x+1}
     P_j_dxyz_i(m,k,1,q,p) = P_new(m,k,q,p) * dble(nx)
    enddo
   enddo
  enddo
 enddo
 ! 
 ! Second contribution is -2 * alpha_i * (x - A_x)^{a_x+1}
 n_new = 1
 call give_poly_i_plus_n_j_xyz(i,j,P_new,iorder_p,n_new)
 do p = 1, ao_prim_num(j)
  do q = 1, ao_prim_num(i)
   do k = 1, 3
    iorder_j_dxyz_i(k,2,q,p) = iorder_p(k,q,p)
    do m = 0, iorder_p(k,q,p)
     ! This corresponds to (x-A_x)^{a_x+1}
     P_j_dxyz_i(m,k,2,q,p) = -2.d0 * ao_expo_ordered_transp(q,i) * P_new(m,k,q,p) 
    enddo
   enddo
  enddo
 enddo


end

subroutine give_poly_i_plus_n_j_xyz(i,j,P_new,iorder_p,n_new)
 implicit none
 include 'utils/constants.include.F'                                                                                                                                  
 integer, intent(in) :: i,j,n_new(3)
 double precision, intent(out) :: P_new(0:max_dim,3,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 integer, intent(out)          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
! BEGIN_DOC
! You change the order of the polynom of phi_i(r) by 
!
! (x-A_x)^{a_x + n_new(1)} 
!
! (y-A_y)^{a_y + n_new(2)}
!
! (z-A_z)^{a_z + n_new(3)}
!
! and then you obtain all the new polynoms for the couple of primitives of the new phi_i * phi_j
!
! 
! END_DOC
 double precision :: P_center(3),fact_p,p_exp
 integer :: dim1,k,p,q,num_i,num_j
 integer :: I_power(3), J_power(3)
 double precision :: I_center(3), J_center(3)
 dim1 = n_pt_max_integrals
 num_i = ao_nucl(i)
 num_j = ao_nucl(j)
 do k = 1, 3
   I_power(k) = ao_power(i,k)
   J_power(k) = ao_power(j,k)
   I_center(k) = nucl_coord(num_i,k)
   J_center(k) = nucl_coord(num_j,k)
 enddo
 do k = 1, 3
  I_power(k) += n_new(k)
  if(I_power(k).lt.0)then
   P_new(:,k,:,:) = 0.d0
   iorder_p(k,:,:) = -1
  endif
 enddo
 do p = 1, ao_prim_num(j)
  do q = 1, ao_prim_num(i)
   do k = 1, 3
    call give_explicit_poly_and_gaussian_x(P_new(0,k,q,p),P_center,p_exp,fact_p,iorder_p(k,q,p),ao_expo_ordered_transp(q,i),ao_expo_ordered_transp(p,j),I_power(k),J_power(k),I_center(k),J_center(k),dim1)
   enddo
  enddo
 enddo
end
