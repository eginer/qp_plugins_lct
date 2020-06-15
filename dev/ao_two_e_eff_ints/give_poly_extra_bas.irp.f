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
