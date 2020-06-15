program ao_two_e_eff_ints
 implicit none
! call test_gauss_ints_aos
! call test_extra_basis
 call test_extra
end
subroutine test_fits
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  integer :: i,nx
  double precision :: x,dx,xmax,fit_1_erf_x,fit_1_erf_x_2
  double precision :: eff_pot_gauss,eff_pot_fit_gauss
  xmax = 10.d0
  nx = 1000
  dx = xmax/dble(nx)
  x = dx
  do i = 1, nx
   write(33,'(100(F16.10,X))')x,(1.d0 - derf(mu_erf*x)),fit_1_erf_x(x),(1.d0 - derf(mu_erf*x))**2.d0,fit_1_erf_x_2(x),eff_pot_gauss(x,mu_erf),eff_pot_fit_gauss(x)
   x += dx
  enddo
end

subroutine test_gauss_ints
 implicit none 
 double precision :: overlap_gauss_r12,weight,r(3),NAI_pol_mult_erf,ERI_erf
 double precision :: D_center(3), C_center(3), A_center(3), B_center(3),alpha,beta,delta,gama,coef
 double precision :: primitive_value_explicit,gauss_a,gauss_b,int_r2
 integer :: power_A(3), power_B(3), power_C(3), power_D(3)
 integer :: ipoint,i,j,n_pt_in
 double precision :: int_erf_num,int_erf_expl
 n_pt_in = n_pt_max_integrals                                                                                            

 power_A = 0
 power_B = 0
 power_C = 0
 power_D = 0
 power_A(1) = 1
 power_B(1) = 1
 power_C(1) = 2
 power_D(1) = 2
 alpha = 1.d0
 beta  = 1.5d0
 gama  = 1.d0
 delta = 4.d0
 A_center = 0.d0
 B_center = 0.d0
 C_center = 0.d0
 D_center = 0.d0

 int_erf_num  = 0.d0
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  gauss_a = primitive_value_explicit(power_A,A_center,alpha,r) ! r1
  gauss_b = primitive_value_explicit(power_B,B_center,beta ,r) ! r1
  int_r2  = NAI_pol_mult_erf(C_center,D_center,power_C,power_D,gama,delta,r,n_pt_in,mu_erf) ! int over r2
  int_erf_num += weight * gauss_a * gauss_b * int_r2 
 enddo
 integer :: a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z
 a_x = power_A(1)
 b_x = power_B(1)
 c_x = power_C(1)
 d_x = power_D(1)

 a_y = power_A(2)
 b_y = power_B(2)
 c_y = power_C(2)
 d_y = power_D(2)

 a_z = power_A(3)
 b_z = power_B(3)
 c_z = power_C(3)
 d_z = power_D(3)

 int_erf_expl = ERI_erf(alpha,beta,delta,gama,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)
 print*,'int_erf_expl = ',int_erf_expl
 print*,'int_erf_num  = ',int_erf_num

end


subroutine test_gauss_ints_bis
 implicit none 
 double precision :: overlap_gauss_r12,weight,r(3),NAI_pol_mult_erf,ERI_erf,weightj,r2(3)
 double precision :: D_center(3), C_center(3), A_center(3), B_center(3),alpha,beta,delta,gama,coef
 double precision :: primitive_value_explicit,gauss_a,gauss_b,int_r2,r_12
 integer :: power_A(3), power_B(3), power_C(3), power_D(3)
 integer :: ipoint,i,j,n_pt_in,jpoint
 double precision :: int_erf_num,int_erf_expl,int_gauss_num,alpha_r12,int_gauss_num_2
 double precision :: general_primitive_integral_gauss
 include 'utils/constants.include.F'

 n_pt_in = n_pt_max_integrals                                                                                            

 power_A = 0
 power_A(1) = 2
 alpha = 1.d0
 A_center = 0.d0
 A_center(1) = 0.5d0

 power_D = 0
 power_C = 0
 power_C(1) = 2
 gama  = 1.d0
 C_center = 0.d0

 int_erf_num  = 0.d0
 int_gauss_num= 0.d0
 int_gauss_num_2 = 0.d0
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  gauss_a = primitive_value_explicit(power_A,A_center,alpha,r) ! r1
  int_r2  = NAI_pol_mult_erf(C_center,C_center,power_C,power_D,0.5d0*gama,0.5d0*gama,r,n_pt_in,mu_erf) ! int over r2
  int_erf_num += weight * gauss_a * int_r2 
  int_r2 = 0.d0
  do i = 1,n_gauss_eff_pot
   alpha_r12 = expo_gauss_eff_pot(i)
   coef      = coef_gauss_eff_pot(i)
   int_r2    = overlap_gauss_r12(r,alpha_r12,C_center,C_center,power_C,power_D,0.5d0*gama,0.5d0*gama)
   int_gauss_num += weight * gauss_a * int_r2 * coef 
  enddo
  int_r2 = 0.d0
  do jpoint = 1, n_points_final_grid
   r2(1) = final_grid_points(1,jpoint)
   r2(2) = final_grid_points(2,jpoint)
   r2(3) = final_grid_points(3,jpoint)
   weightj = final_weight_at_r_vector(jpoint)
   r_12 = (r(1) - r2(1))**2 + (r(2) - r2(2))**2 + (r(3) - r2(3))**2 
   gauss_b = primitive_value_explicit(power_C,C_center,gama,r2) ! r1
   do i = 1,n_gauss_eff_pot
    alpha_r12 = expo_gauss_eff_pot(i)
    if(alpha_r12 * r_12.gt.20.d0)cycle
    coef      = coef_gauss_eff_pot(i)
    int_r2   += gauss_b * dexp(-alpha_r12*r_12) * coef * weightj
   enddo
  enddo
  int_gauss_num_2 += weight * gauss_a * int_r2 
 enddo
 integer :: dim1
 dim1 = n_pt_max_integrals
 double precision, allocatable :: P_new(:,:),Q_new(:,:)
 integer :: iorder_p(3), iorder_q(3)
 double precision :: P_center(3),fact_p,pp,p_inv
 double precision :: Q_center(3),fact_q,qq,q_inv
 double precision :: erf_int,gauss_int
 allocate(P_new(0:max_dim,3),Q_new(0:max_dim,3))

 P_new = 0.d0
 P_new(power_A(1),1)=1.d0
 P_new(power_A(2),2)=1.d0
 P_new(power_A(3),3)=1.d0
 P_center = A_center
 pp = alpha
 p_inv = 1.d0/pp
 iorder_p = power_A
 fact_p = 1.d0

 Q_new = 0.d0
 Q_new(power_C(1),1)=1.d0
 Q_new(power_C(2),2)=1.d0
 Q_new(power_C(3),3)=1.d0
 Q_center = C_center
 qq = gama
 q_inv = 1.d0/qq
 iorder_q = power_C
 fact_q = 1.d0

 gauss_int =  general_primitive_integral_gauss(dim1,                 &
      P_new,P_center,fact_p,pp,p_inv,iorder_p,                   &
      Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
 print*,'int_erf_num  = ',int_erf_num
! print*,'erf_int      = ',erf_int

 print*,'int_gauss_num= ',int_gauss_num
 print*,'int_gauss_num_2',int_gauss_num_2
 print*,'gauss_int    = ',gauss_int

end


