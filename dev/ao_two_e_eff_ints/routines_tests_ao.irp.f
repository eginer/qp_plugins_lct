subroutine test_gauss_ints_ao
 implicit none 
 double precision :: overlap_gauss_r12,weight,r(3),weightj,r2(3)
 double precision :: D_center(3), C_center(3), A_center(3), B_center(3),alpha,beta,delta,gama,coef
 double precision :: primitive_value_explicit,gauss_a,gauss_b,int_r2,r_12
 integer :: power_A(3), power_B(3), power_C(3), power_D(3)
 integer :: ipoint,i,j,n_pt_in,jpoint
 double precision :: int_gauss_num,alpha_r12,gauss_b
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

 int_gauss_num= 0.d0
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  gauss_a = primitive_value_explicit(power_A,A_center,alpha,r) ! r1
  int_r2 = 0.d0
  do i = 1,n_gauss_eff_pot
   alpha_r12 = expo_gauss_eff_pot(i)
   coef      = coef_gauss_eff_pot(i)
   int_r2    = overlap_gauss_r12(r,alpha_r12,C_center,C_center,power_C,power_D,0.5d0*gama,0.5d0*gama)
   int_gauss_num += weight * gauss_a * int_r2 * coef 
  enddo
 enddo
 integer :: dim1
 dim1 = n_pt_max_integrals
 double precision, allocatable :: P_new(:,:),Q_new(:,:)
 integer :: iorder_p(3), iorder_q(3)
 double precision :: P_center(3),fact_p,pp,p_inv
 double precision :: Q_center(3),fact_q,qq,q_inv
 double precision :: gauss_int
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

 print*,'int_gauss_num= ',int_gauss_num
 print*,'gauss_int    = ',gauss_int

end


