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
   if(isnan(ao_prod_in_r))then
    print*,'isnan(ao_prod_in_r)',isnan(ao_prod_in_r)
    print*,'p,q',p,q
    print*,r
    print*,P_center(:,q,p)
    print*,dr,p_exp(q,p)
    write(*,'(100(F16.10,X))')pol_tot,gauss_r,fact_p(q,p),coef_prod(q,p)
    pause
   endif
  enddo
 enddo
end

double precision function ao_prod_in_r_bis(r,prim_num_i,prim_num_j,P_x,P_y,P_z,P_center,p_exp,fact_p,iorder_p,coef_prod)
  include 'utils/constants.include.F'
 double precision, intent(in) :: r(3)
 integer, intent(in) :: prim_num_i,prim_num_j
 double precision, intent(in) :: P_x(0:max_dim,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 double precision, intent(in) :: P_y(0:max_dim,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 double precision, intent(in) :: P_z(0:max_dim,ao_prim_num_max,ao_prim_num_max) ! new polynom for each couple of prim
 double precision, intent(in) :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision, intent(in) :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision, intent(in) :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer, intent(in)          :: iorder_p(3,ao_prim_num_max,ao_prim_num_max) ! order of the polynoms for each couple of prim
 double precision, intent(in) :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision :: pol(3),dx(3),dr,gauss_r,pol_tot
 integer :: m,p,q,n
 ao_prod_in_r_bis = 0.d0
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
   do n = 0, iorder_p(1,q,p)
    pol(1) += dx(1)**dble(n) * P_x(n,q,p)
   enddo
   do n = 0, iorder_p(2,q,p)
    pol(2) += dx(2)**dble(n) * P_y(n,q,p)
   enddo
   do n = 0, iorder_p(3,q,p)
    pol(3) += dx(3)**dble(n) * P_z(n,q,p)
   enddo
   pol_tot = pol(1) * pol(2) * pol(3)
   ao_prod_in_r_bis += pol_tot * gauss_r * fact_p(q,p) * coef_prod(q,p)
!   write(*,'(100(F16.10,X))')pol_tot,gauss_r,fact_p(q,p),coef_prod(q,p)
   if(isnan(ao_prod_in_r_bis))then
    print*,'isnan(ao_prod_in_r_bis)',isnan(ao_prod_in_r_bis)
    print*,'p,q',p,q
    print*,r
    print*,P_center(:,q,p)
    print*,dr,p_exp(q,p)
    write(*,'(100(F16.10,X))')pol_tot,gauss_r,fact_p(q,p),coef_prod(q,p)
    pause
   endif
  enddo
 enddo
! print*,"ao_prod_in_r_bis = ",ao_prod_in_r_bis
end

subroutine ao_xyz_i_prod_j_in_r(r,i,j,aos_xyzi_j)
  include 'utils/constants.include.F'
  implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: i,j
 double precision, intent(out) :: aos_xyzi_j(3)
 double precision :: P_new0(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_new1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer          :: iorder_p0(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer          :: iorder_p1(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 double precision :: coef_prod(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision :: ao_prod_in_r
 integer :: prim_num_j,prim_num_i,k
 aos_xyzi_j = 0.d0
 call give_poly_xyzi_j(i,j,P_new0,P_new1,P_center,p_exp,fact_p,iorder_p0,iorder_p1,coef_prod)
 prim_num_i = ao_prim_num(i)
 prim_num_j = ao_prim_num(j)
 aos_xyzi_j = 0.d0
 double precision :: tmp
 double precision :: aos_array(ao_num),center_ao(3)
 integer :: num_ao
 call give_all_aos_at_r(r,aos_array)
 num_ao = ao_nucl(i)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 do k = 1,3
  tmp = ao_prod_in_r(r,prim_num_i,prim_num_j,P_new0(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p0(1,1,1,k),coef_prod)
  aos_xyzi_j(k) = tmp
  tmp = ao_prod_in_r(r,prim_num_i,prim_num_j,P_new1(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p1(1,1,1,k),coef_prod)
  aos_xyzi_j(k) += tmp
 enddo
end

subroutine ao_dxyz_i_prod_j_in_r(r,i,j,aos_dxyzi_j)
  include 'utils/constants.include.F'
  implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: i,j
 double precision :: aos_dxyzi_j(3)
 double precision :: P_new0(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_new1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer          :: iorder_p0(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer          :: iorder_p1(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 double precision :: coef_prod0(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision :: coef_prod1(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision :: ao_prod_in_r
 integer :: prim_num_j,prim_num_i,k
 aos_dxyzi_j = 0.d0
 call give_poly_dxyzi_j(i,j,P_new0,P_new1,P_center,p_exp,fact_p,iorder_p0,iorder_p1,coef_prod0,coef_prod1)
 prim_num_i = ao_prim_num(i)
 prim_num_j = ao_prim_num(j)
 aos_dxyzi_j = 0.d0
 do k = 1, 3
  aos_dxyzi_j(k)  = ao_prod_in_r(r,prim_num_i,prim_num_j,P_new0(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p0(1,1,1,k),coef_prod0)
  aos_dxyzi_j(k) += ao_prod_in_r(r,prim_num_i,prim_num_j,P_new1(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p1(1,1,1,k),coef_prod1)
 enddo
end

subroutine ao_xyz_dxyz_i_prod_j_in_r(r,i,j,aos_dxyzi_j)
  include 'utils/constants.include.F'
  implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: i,j
 double precision, intent(out) :: aos_dxyzi_j(3)
 double precision :: P_new0(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_new1(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_new2(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_new3(0:max_dim,3,ao_prim_num_max,ao_prim_num_max,3) ! new polynom for each couple of prim
 double precision :: P_center(3,ao_prim_num_max,ao_prim_num_max) ! new center for each couple of prim
 double precision :: p_exp(ao_prim_num_max,ao_prim_num_max) ! new gaussian exponents for each couple of prim
 double precision :: fact_p(ao_prim_num_max,ao_prim_num_max) ! factor for each couple of primitive 
 integer          :: iorder_p0(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer          :: iorder_p1(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer          :: iorder_p2(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 integer          :: iorder_p3(3,ao_prim_num_max,ao_prim_num_max,3) ! order of the polynoms for each couple of prim
 double precision :: coef_prod0(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision :: coef_prod1(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision :: coef_prod2(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 
 double precision :: coef_prod3(ao_prim_num_max,ao_prim_num_max) ! produc of coef for each couple of primitive 

 double precision :: ao_prod_in_r
 integer :: prim_num_j,prim_num_i,k
 aos_dxyzi_j = 0.d0
 call give_poly_xyz_dxyzi_j(i,j,P_new0,P_new1,P_new2,P_new3,P_center,p_exp,fact_p,iorder_p0,iorder_p1,iorder_p2,iorder_p3,coef_prod0,coef_prod1, coef_prod2, coef_prod3)
 prim_num_i = ao_prim_num(i)
 prim_num_j = ao_prim_num(j)
 aos_dxyzi_j = 0.d0
 do k = 1, 3
  aos_dxyzi_j(k)  = ao_prod_in_r(r,prim_num_i,prim_num_j,P_new0(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p0(1,1,1,k),coef_prod0)
  aos_dxyzi_j(k) += ao_prod_in_r(r,prim_num_i,prim_num_j,P_new1(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p1(1,1,1,k),coef_prod1)
  aos_dxyzi_j(k) += ao_prod_in_r(r,prim_num_i,prim_num_j,P_new2(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p2(1,1,1,k),coef_prod2)
  aos_dxyzi_j(k) += ao_prod_in_r(r,prim_num_i,prim_num_j,P_new3(0,1,1,1,k),P_center,p_exp,fact_p,iorder_p3(1,1,1,k),coef_prod3)
 enddo
end



double precision function phi_ao_plus_n(r,i_ao,n_new,ixyz)
 implicit none
 integer, intent(in) :: i_ao,n_new,ixyz
 double precision, intent(in) :: r(3)
 double precision :: center_ao(3),beta
 double precision :: accu,dr(3),r2,pol_usual(3)
 integer :: m,power_ao(3),num_ao
 power_ao(1:3)= ao_power(i_ao,1:3) 

 power_ao(ixyz) += n_new

 num_ao = ao_nucl(i_ao)
 center_ao(1:3) = nucl_coord(num_ao,1:3)
 dr(1) = (r(1) - center_ao(1))
 dr(2) = (r(2) - center_ao(2))
 dr(3) = (r(3) - center_ao(3))
 r2 = 0.d0
 do m = 1, 3
  r2 += dr(m)*dr(m)
 enddo
 dr(1) = dr(1)**dble(power_ao(1))
 dr(2) = dr(2)**dble(power_ao(2))
 dr(3) = dr(3)**dble(power_ao(3))
 ! computes the gaussian part 
 accu = 0.d0
 do m=1,ao_prim_num(i_ao)
   beta = ao_expo_ordered_transp(m,i_ao)
   if(dabs(beta*r2).gt.50.d0)cycle
   accu += ao_coef_normalized_ordered_transp(m,i_ao) * dexp(-beta*r2)
 enddo
 ! computes the polynom part
 phi_ao_plus_n = accu * dr(1) * dr(2) * dr(3)
end
