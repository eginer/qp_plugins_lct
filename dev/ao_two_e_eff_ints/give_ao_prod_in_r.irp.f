double precision function ao_prod_in_r(r,prim_num_i,prim_num_j,P_new,P_center,p_exp,fact_p,iorder_p,coef_prod)
  include 'utils/constants.include.F'
 double precision, intent(in) :: r(3)
 integer, intent(in) :: prim_num_i,prim_num_j
 double precision, intent(in) :: P_new(0:max_dim,3,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 double precision, intent(in) :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(in) :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(in) :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer, intent(in)          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision, intent(in) :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
! call give_poly_ij(i,j,P_new,P_center,p_exp,fact_p,iorder_p)
 double precision :: pol(3),dx(3),dr,gauss_r,pol_tot
 integer :: m,p,q,n
 ao_prod_in_r = 0.d0
 do p = 1, prim_num_j
  do q = 1, prim_num_i
   pol = 0.d0
   dr = 0.d0
   do m = 1, 3
    dx(m) = r(m) - P_center(m,q,p)
    dr += dx(m) * dx(m)
   enddo
   if(dabs(dr * p_exp(q,p)).gt.50.d0)then
    gauss_r = 0.d0
   else
    gauss_r = dexp(-dr * p_exp(q,p))
   endif
   do m = 1, 3
    do n = 0, iorder_p(m,q,p)
     pol(m) += dx(m)**dble(n) * P_new(n,m,q,p)
    enddo
   enddo
   pol_tot = pol(1) * pol(2) * pol(3)
   ao_prod_in_r += pol_tot * gauss_r * fact_p(q,p) * coef_prod(q,p)
  enddo
 enddo
end
