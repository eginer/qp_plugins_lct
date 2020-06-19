subroutine give_poly_ij(i,j,P_new,P_center,p_exp,fact_p,iorder_p,coef_prod)
 implicit none
  include 'utils/constants.include.F'
 BEGIN_DOC
! provide all the polynoms, order of polynoms, exponential factors, centers for all couple of primitive of two AO i and j
!
! fact_k(m,n) = exponential factor for the "mth" primitive of AO i and "nth" primitive of AO j
!
! P_new(0:max_dim,3,m,n) = polynom for the "mth" primitive of AO i and "nth" primitive of AO j
!
! p_exp(m,n)  = gaussian exponent   for the "mth" primitive of AO i and "nth" primitive of AO j
!
! iorder(3,m,n) = order of the polynomx for the "mth" primitive of AO i and "nth" primitive of AO j
!
! coef_prod(m,n) = product of coefficients for the "mth" primitive of AO i and "nth" primitive of AO j
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(out) :: P_new(0:max_dim,3,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 double precision, intent(out) :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(out) :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(out) :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer, intent(out)          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 integer :: p,q,I_power(3),J_power(3),num_i,num_j,dim1
 double precision :: I_center(3),J_center(3),coef1,coef2
 dim1 = n_pt_max_integrals

 num_i = ao_nucl(i)
 num_j = ao_nucl(j)

 do p = 1, 3
   I_power(p) = ao_power(i,p)
   J_power(p) = ao_power(j,p)
   I_center(p) = nucl_coord(num_i,p)
   J_center(p) = nucl_coord(num_j,p)
 enddo
 do p = 1, ao_prim_num(j)
   coef1 = ao_coef_normalized_ordered_transp(p,j)
   do q = 1, ao_prim_num(i)
    coef2 = coef1*ao_coef_normalized_ordered_transp(q,i)
    coef_prod(q,p) = coef2 
    call give_explicit_poly_and_gaussian(P_new(0,1,q,p),P_center(1,q,p),p_exp(q,p),fact_p(q,p),iorder_p(1,q,p),&
        ao_expo_ordered_transp(q,i),ao_expo_ordered_transp(p,j),                 &
        I_power,J_power,I_center,J_center,dim1)
   enddo
 enddo

end

subroutine give_poly_i_plus_n_j(i,j,P_new,P_center,p_exp,fact_p,iorder_p,coef_prod,n_new,ixyz)
 implicit none
  include 'utils/constants.include.F'
 BEGIN_DOC
! provide all the polynoms, order of polynoms, exponential factors, centers for all couple of primitive of AO i and j
!
! fact_k(m,n) = exponential factor for the "mth" primitive of AO i and "nth" primitive of AO j
!
! P_new(0:max_dim,3,m,n) = polynom for the "mth" primitive of AO i and "nth" primitive of AO j
!
! p_exp(m,n)  = gaussian exponent   for the "mth" primitive of AO i and "nth" primitive of AO j
!
! iorder(3,m,n) = order of the polynomx for the "mth" primitive of AO i and "nth" primitive of AO j
!
! coef_prod(m,n) = product of coefficients for the "mth" primitive of AO i and "nth" primitive of AO j
 END_DOC
 integer, intent(in) :: i,j,n_new,ixyz
 double precision, intent(out) :: P_new(0:max_dim,3,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 double precision, intent(out) :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(out) :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(out) :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer, intent(out)          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 integer :: p,q,I_power(3),J_power(3),num_i,num_j,dim1
 double precision :: I_center(3),J_center(3),coef1,coef2
 dim1 = n_pt_max_integrals

 num_i = ao_nucl(i)
 num_j = ao_nucl(j)

 do p = 1, 3
   I_power(p) = ao_power(i,p)
   J_power(p) = ao_power(j,p)
   I_center(p) = nucl_coord(num_i,p)
   J_center(p) = nucl_coord(num_j,p)
 enddo
 I_power(ixyz) += n_new
 if(I_power(ixyz).lt.0)then
  P_new = 0.d0
  P_center = 0.d0
  p_exp = 1.d10
  fact_p = 0.d0
  iorder_p = -1
  coef_prod = 0.d0
  return
 endif
 do p = 1, ao_prim_num(j)
   coef1 = ao_coef_normalized_ordered_transp(p,j)
   do q = 1, ao_prim_num(i)
    coef2 = coef1*ao_coef_normalized_ordered_transp(q,i)
    coef_prod(q,p) = coef2 
    call give_explicit_poly_and_gaussian(P_new(0,1,q,p),P_center(1,q,p),p_exp(q,p),fact_p(q,p),iorder_p(1,q,p),&
        ao_expo_ordered_transp(q,i),ao_expo_ordered_transp(p,j),                 &
        I_power,J_power,I_center,J_center,dim1)
   enddo
 enddo

end

subroutine give_poly_xyzi_j(i,j,P_new0,P_new1,P_center,p_exp,fact_p,iorder_p0,iorder_p1,coef_prod)
 implicit none
  include 'utils/constants.include.F'
 BEGIN_DOC
! provide all the polynoms, order of polynoms, exponential factors, centers for all couple of primitive for {x,y,z}* AO i and AO j
!
! fact_k(m,n) = exponential factor for the "mth" primitive of AO i and "nth" primitive of AO j
!
! P_new(0:max_dim,3,m,n) = polynom for the "mth" primitive of AO i and "nth" primitive of AO j
!
! p_exp(m,n)  = gaussian exponent   for the "mth" primitive of AO i and "nth" primitive of AO j
!
! iorder(3,m,n) = order of the polynomx for the "mth" primitive of AO i and "nth" primitive of AO j
!
! coef_prod(m,n) = product of coefficients for the "mth" primitive of AO i and "nth" primitive of AO j
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(out) :: P_new0(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_new1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(out) :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(out) :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer, intent(out)          :: iorder_p0(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer, intent(out)          :: iorder_p1(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision :: P_new_ij(0:max_dim,3,ao_prim_num_max,ao_prim_num_max)
 double precision :: P_center_ij(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: p_exp_ij(ao_prim_num_max,ao_prim_num_max) 
 double precision :: fact_p_ij(ao_prim_num_max,ao_prim_num_max) 
 integer          :: iorder_p_ij(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: coef_prod_ij(ao_prim_num_max,ao_prim_num_max) 

 double precision :: P_new_ij_p1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max)
 double precision :: P_center_ij_p1(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: p_exp_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 double precision :: fact_p_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 integer          :: iorder_p_ij_p1(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: coef_prod_ij_p1(ao_prim_num_max,ao_prim_num_max) 

 P_new0 = 0.d0
 P_new1 = 0.d0
 call give_poly_ij(i,j,P_new_ij,P_center_ij,p_exp_ij,fact_p_ij,iorder_p_ij,coef_prod_ij)
 P_center = P_center_ij
 fact_p   = fact_p_ij
 coef_prod= coef_prod_ij
 p_exp    = p_exp_ij

 integer :: p,q,k,m,n,num_ao
 double precision :: center_ao(1:3)
 integer :: n_new,ixyz
 num_ao = ao_nucl(i)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 ! initialize the polynoms to the polynoms of phi_i * phi_j
 do k = 1, 3
  iorder_p0(:,:,:,k) = iorder_p_ij(:,:,:)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     do m = 0, iorder_p_ij(n,q,p)
      P_new0(m,n,q,p,k) = P_new_ij(m,n,q,p) 
     enddo
    enddo
   enddo
  enddo
 enddo
 ! Then, x * (x - A_x)^a_x = A_x (x - A_x)^a_x + 1 * (x - A_x)^{a_x+1}
 ! So you have to multiply by A_x the polynom on x of phi_i * phi_j 
 do k = 1, 3
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    n = k
    do m = 0, iorder_p_ij(k,q,p)
     P_new0(m,n,q,p,k) = P_new_ij(m,k,q,p) * center_ao(k)
    enddo
   enddo
  enddo
 enddo
 ! Then you have to add the new polynoms with (x-A_x)^{a_x+1}
 n_new = 1
 do k = 1, 3
  ixyz = k
  call give_poly_i_plus_n_j(i,j,P_new_ij_p1,P_center_ij_p1,p_exp_ij_p1,fact_p_ij_p1,iorder_p_ij_p1,coef_prod_ij_p1,n_new,ixyz)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     iorder_p1(n,q,p,k) = iorder_p_ij_p1(n,q,p)
     do m = 0, iorder_p_ij_p1(n,q,p)
      P_new1(m,n,q,p,k) = P_new_ij_p1(m,n,q,p) 
     enddo
    enddo
   enddo
  enddo
 enddo

end

subroutine give_poly_dxyzi_j(i,j,P_new0,P_new1,P_center,p_exp,fact_p,iorder_p0,iorder_p1,coef_prod0,coef_prod1)
 implicit none
  include 'utils/constants.include.F'
 BEGIN_DOC
! provide all the polynoms, order of polynoms, exponential factors, centers for all couple of primitive for {x,y,z}* AO i and AO j
!
! fact_k(m,n) = exponential factor for the "mth" primitive of AO i and "nth" primitive of AO j
!
! P_new(0:max_dim,3,m,n) = polynom for the "mth" primitive of AO i and "nth" primitive of AO j
!
! p_exp(m,n)  = gaussian exponent   for the "mth" primitive of AO i and "nth" primitive of AO j
!
! iorder(3,m,n) = order of the polynomx for the "mth" primitive of AO i and "nth" primitive of AO j
!
! coef_prod(m,n) = product of coefficients for the "mth" primitive of AO i and "nth" primitive of AO j
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(out) :: P_new0(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_new1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(out) :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(out) :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer, intent(out)          :: iorder_p0(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer, intent(out)          :: iorder_p1(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: coef_prod0(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision, intent(out) :: coef_prod1(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision :: P_new_ij(0:max_dim,3,ao_prim_num_max,ao_prim_num_max)
 double precision :: P_center_ij(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: p_exp_ij(ao_prim_num_max,ao_prim_num_max) 
 double precision :: fact_p_ij(ao_prim_num_max,ao_prim_num_max) 
 integer          :: iorder_p_ij(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: coef_prod_ij(ao_prim_num_max,ao_prim_num_max) 

 double precision :: P_new_ij_p1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max)
 double precision :: P_center_ij_p1(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: p_exp_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 double precision :: fact_p_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 integer          :: iorder_p_ij_p1(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: coef_prod_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 
 double precision :: coef1,coef2
 integer :: nx,p,q,k,n_new,ixyz,n,m
 do p = 1, ao_prim_num(j)
   coef1 = ao_coef_normalized_ordered_transp(p,j)
   do q = 1, ao_prim_num(i)
    coef2 = coef1*ao_coef_normalized_ordered_transp(q,i)
    coef_prod0(q,p) = coef2 
    coef_prod1(q,p) = -2.d0 * coef2 * ao_expo_ordered_transp(q,i)
   enddo
 enddo

 call give_poly_ij(i,j,P_new_ij,P_center_ij,p_exp_ij,fact_p_ij,iorder_p_ij,coef_prod_ij_p1)
 P_center = P_center_ij
 fact_p   = fact_p_ij
 p_exp    = p_exp_ij

 P_new0 = 0.d0
 P_new1 = 0.d0
 n_new = -1
 do k = 1, 3
  nx = ao_power(i,k)
  ixyz = k
  call give_poly_i_plus_n_j(i,j,P_new_ij_p1,P_center_ij_p1,p_exp_ij_p1,fact_p_ij_p1,iorder_p_ij_p1,coef_prod_ij_p1,n_new,ixyz)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     iorder_p0(n,q,p,k) = iorder_p_ij_p1(n,q,p)
     do m = 0, iorder_p_ij_p1(n,q,p)
      P_new0(m,n,q,p,k) = P_new_ij_p1(m,n,q,p) * dble(nx)
     enddo
    enddo
   enddo
  enddo
 enddo

 n_new = 1
 do k = 1, 3
  nx = ao_power(i,k)
  ixyz = k
  call give_poly_i_plus_n_j(i,j,P_new_ij_p1,P_center_ij_p1,p_exp_ij_p1,fact_p_ij_p1,iorder_p_ij_p1,coef_prod_ij_p1,n_new,ixyz)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     iorder_p1(n,q,p,k) = iorder_p_ij_p1(n,q,p)
     do m = 0, iorder_p_ij_p1(n,q,p)
      P_new1(m,n,q,p,k) = P_new_ij_p1(m,n,q,p) 
     enddo
    enddo
   enddo
  enddo
 enddo

end

subroutine give_poly_xyz_dxyzi_j(i,j,P_new0,P_new1,P_new2,P_new3,P_center,p_exp,fact_p,iorder_p0,iorder_p1,iorder_p2,iorder_p3,coef_prod0,coef_prod1, coef_prod2, coef_prod3)
 implicit none
  include 'utils/constants.include.F'
 BEGIN_DOC
! provide all the polynoms, order of polynoms, exponential factors, centers for all couple of primitive for {x,y,z}* AO i and AO j
!
! fact_k(m,n) = exponential factor for the "mth" primitive of AO i and "nth" primitive of AO j
!
! P_new(0:max_dim,3,m,n) = polynom for the "mth" primitive of AO i and "nth" primitive of AO j
!
! p_exp(m,n)  = gaussian exponent   for the "mth" primitive of AO i and "nth" primitive of AO j
!
! iorder(3,m,n) = order of the polynomx for the "mth" primitive of AO i and "nth" primitive of AO j
!
! coef_prod(m,n) = product of coefficients for the "mth" primitive of AO i and "nth" primitive of AO j
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(out) :: P_new0(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_new1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_new2(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_new3(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision, intent(out) :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(out) :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(out) :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer, intent(out)          :: iorder_p0(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer, intent(out)          :: iorder_p1(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer, intent(out)          :: iorder_p2(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer, intent(out)          :: iorder_p3(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 double precision, intent(out) :: coef_prod0(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision, intent(out) :: coef_prod1(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision, intent(out) :: coef_prod2(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision, intent(out) :: coef_prod3(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision :: P_new_ij(0:max_dim,3,ao_prim_num_max,ao_prim_num_max)
 double precision :: P_center_ij(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: p_exp_ij(ao_prim_num_max,ao_prim_num_max) 
 double precision :: fact_p_ij(ao_prim_num_max,ao_prim_num_max) 
 integer          :: iorder_p_ij(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: coef_prod_ij(ao_prim_num_max,ao_prim_num_max) 

 double precision :: P_new_ij_p1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max)
 double precision :: P_center_ij_p1(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: p_exp_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 double precision :: fact_p_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 integer          :: iorder_p_ij_p1(3,ao_prim_num_max,ao_prim_num_max) 
 double precision :: coef_prod_ij_p1(ao_prim_num_max,ao_prim_num_max) 
 
 double precision :: coef1,coef2,center_ao(3)
 integer :: nx,p,q,k,n_new,ixyz,n,m,num_ao
 P_new0 = 0.d0
 P_new1 = 0.d0
 P_new2 = 0.d0
 P_new3 = 0.d0
 num_ao = ao_nucl(i)
 center_ao(1:3) = nucl_coord(num_ao,1:3)


 
 ! a_x * Ax * (x - Ax)^{a_x-1}
 n_new = -1
 do k = 1, 3
  nx = ao_power(i,k)
  ixyz = k
  call give_poly_i_plus_n_j(i,j,P_new_ij_p1,P_center_ij_p1,p_exp_ij_p1,fact_p_ij_p1,iorder_p_ij_p1,coef_prod_ij_p1,n_new,ixyz)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     iorder_p0(n,q,p,k) = iorder_p_ij_p1(n,q,p)
     do m = 0, iorder_p_ij_p1(n,q,p)
      P_new0(m,n,q,p,k) = P_new_ij_p1(m,n,q,p) 
     enddo
    enddo
   enddo
  enddo
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
     n = k 
     do m = 0, iorder_p_ij_p1(n,q,p)
      P_new0(m,n,q,p,k) = P_new_ij_p1(m,n,q,p) * dble(nx) * center_ao(k)
     enddo
   enddo
  enddo
 enddo

 ! a_x * (x - Ax)^a_x
 call give_poly_ij(i,j,P_new_ij,P_center_ij,p_exp_ij,fact_p_ij,iorder_p_ij,coef_prod_ij_p1)
 P_center = P_center_ij
 fact_p   = fact_p_ij
 p_exp    = p_exp_ij
 coef_prod1 = coef_prod_ij_p1 ! polynom of order ax
 coef_prod0 = coef_prod_ij_p1 ! polynom of order ax-1

 ! initialize the polynom of order ax to the polynoms of phi_i * phi_j
 do k = 1, 3
  iorder_p1(:,:,:,k) = iorder_p_ij(:,:,:)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     do m = 0, iorder_p_ij(n,q,p)
      P_new1(m,n,q,p,k) = P_new_ij(m,n,q,p) 
     enddo
    enddo
   enddo
  enddo
 enddo
 ! So you have to multiply by ax the polynom on x of phi_i * phi_j 
 do k = 1, 3
  nx = ao_power(i,k)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    n = k
    do m = 0, iorder_p_ij(k,q,p)
     P_new1(m,n,q,p,k) = P_new_ij(m,k,q,p) * dble(nx)
    enddo
   enddo
  enddo
 enddo

 ! - 2 alpha_i * Ax * (x - Ax)^{ax+1}
 do p = 1, ao_prim_num(j)
   coef1 = ao_coef_normalized_ordered_transp(p,j)
   do q = 1, ao_prim_num(i)
    coef2 = coef1*ao_coef_normalized_ordered_transp(q,i)
    coef_prod2(q,p) = -2.d0 * coef2 * ao_expo_ordered_transp(q,i)
   enddo
 enddo

 n_new = 1
 do k = 1, 3
  nx = ao_power(i,k)
  ixyz = k
  call give_poly_i_plus_n_j(i,j,P_new_ij_p1,P_center_ij_p1,p_exp_ij_p1,fact_p_ij_p1,iorder_p_ij_p1,coef_prod_ij_p1,n_new,ixyz)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     iorder_p2(n,q,p,k) = iorder_p_ij_p1(n,q,p)
     do m = 0, iorder_p_ij_p1(n,q,p)
      P_new2(m,n,q,p,k) = P_new_ij_p1(m,n,q,p)
     enddo
    enddo
   enddo
  enddo
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    n = k
    do m = 0, iorder_p_ij_p1(n,q,p)
     P_new2(m,n,q,p,k) = P_new_ij_p1(m,n,q,p) * center_ao(k)
    enddo
   enddo
  enddo
 enddo

 ! - 2 alpha_i * (x - Ax)^{ax+2}
 do p = 1, ao_prim_num(j)
   coef1 = ao_coef_normalized_ordered_transp(p,j)
   do q = 1, ao_prim_num(i)
    coef2 = coef1*ao_coef_normalized_ordered_transp(q,i)
    coef_prod3(q,p) = -2.d0 * coef2 * ao_expo_ordered_transp(q,i)
   enddo
 enddo

 n_new = 2
 do k = 1, 3
  ixyz = k
  call give_poly_i_plus_n_j(i,j,P_new_ij_p1,P_center_ij_p1,p_exp_ij_p1,fact_p_ij_p1,iorder_p_ij_p1,coef_prod_ij_p1,n_new,ixyz)
  do p = 1, ao_prim_num(j)
   do q = 1, ao_prim_num(i)
    do n = 1, 3
     iorder_p3(n,q,p,k) = iorder_p_ij_p1(n,q,p)
     do m = 0, iorder_p_ij_p1(n,q,p)
      P_new3(m,n,q,p,k) = P_new_ij_p1(m,n,q,p) 
     enddo
    enddo
   enddo
  enddo
 enddo

end

