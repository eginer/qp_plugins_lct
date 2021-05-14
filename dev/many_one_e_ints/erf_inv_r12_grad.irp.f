subroutine phi_j_erf_mu_r_dxyz_phi(i,j,mu_in, C_center, dxyz_ints)
 implicit none
 BEGIN_DOC
! dxyz_ints(1/2/3) = int dr phi_j(r) [erf(mu  |r - C|)/|r-C|]  d/d(x/y/z) phi_i(r)
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: mu_in, C_center(3)
 double precision, intent(out):: dxyz_ints(3)
 integer :: num_A,power_A(3), num_b, power_B(3),power_B_tmp(3)
 double precision :: alpha, beta, A_center(3), B_center(3),contrib,NAI_pol_mult_erf
 integer :: n_pt_in,l,m,mm
 n_pt_in = n_pt_max_integrals
 ! j 
 num_A = ao_nucl(j)
 power_A(1:3)= ao_power(j,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 ! i 
 num_B = ao_nucl(i)
 power_B(1:3)= ao_power(i,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)

 dxyz_ints = 0.d0
 do l=1,ao_prim_num(j)
  alpha = ao_expo_ordered_transp(l,j)
  do m=1,ao_prim_num(i)
    beta = ao_expo_ordered_transp(m,i)
    do mm = 1, 3
     ! (d/dx phi_i ) * phi_j 
     ! d/dx * (x - B_x)^b_x exp(-beta * (x -B_x)^2)= [b_x * (x - B_x)^(b_x - 1) - 2 beta * (x - B_x)^(b_x + 1)] exp(-beta * (x -B_x)^2)
     !
     ! first contribution :: b_x (x - B_x)^(b_x-1) :: integral with b_x=>b_x-1 multiplied by b_x
     power_B_tmp = power_B
     power_B_tmp(mm) += -1
     contrib = NAI_pol_mult_erf(A_center,B_center,power_A,power_B_tmp,alpha,beta,C_center,n_pt_in,mu_in)  
     dxyz_ints(mm) += contrib * dble(power_B(mm)) * ao_coef_normalized_ordered_transp(l,j)             &
                                                  * ao_coef_normalized_ordered_transp(m,i) 
     ! second contribution ::  - 2 beta * (x - B_x)^(b_x + 1) :: integral with b_x=> b_x+1 multiplied by -2 * beta
     power_B_tmp = power_B
     power_B_tmp(mm) += 1
     contrib = NAI_pol_mult_erf(A_center,B_center,power_A,power_B_tmp,alpha,beta,C_center,n_pt_in,mu_in)  
     dxyz_ints(mm) += contrib * (-2.d0 * beta )  * ao_coef_normalized_ordered_transp(l,j)             &
                                                 * ao_coef_normalized_ordered_transp(m,i) 
    enddo
  enddo
 enddo
end



