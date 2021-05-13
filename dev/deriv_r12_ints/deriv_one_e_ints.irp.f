subroutine phi_j_erf_mu_r_xyz_phi(i,j,mu_in, C_center, xyz_ints)
 implicit none
 BEGIN_DOC
! xyz_ints(1/2/3) = int dr phi_j(r) [erf(mu  |r - C|)/|r-C|]  x/y/z phi_i(r)
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: mu_in, C_center(3)
 double precision, intent(out):: xyz_ints(3)
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

 xyz_ints = 0.d0
 do l=1,ao_prim_num(j)
  alpha = ao_expo_ordered_transp(l,j)
  do m=1,ao_prim_num(i)
    beta = ao_expo_ordered_transp(m,i)
    do mm = 1, 3
     ! (x phi_i ) * phi_j 
     ! x * (x - B_x)^b_x = b_x (x - B_x)^b_x + 1 * (x - B_x)^{b_x+1}
     !
     ! first contribution :: B_x (x - B_x)^b_x :: usual integral multiplied by B_x
     power_B_tmp = power_B
     contrib = NAI_pol_mult_erf(A_center,B_center,power_A,power_B_tmp,alpha,beta,C_center,n_pt_in,mu_in)  
     xyz_ints(mm) += contrib * B_center(mm) * ao_coef_normalized_ordered_transp(l,j)             &
                                            * ao_coef_normalized_ordered_transp(m,i) 
     ! second contribution :: 1 * (x - B_x)^(b_x+1) :: integral with b_x=>b_x+1 
     power_B_tmp(mm) += 1
     contrib = NAI_pol_mult_erf(A_center,B_center,power_A,power_B_tmp,alpha,beta,C_center,n_pt_in,mu_in)  
     xyz_ints(mm) += contrib * 1.d0        * ao_coef_normalized_ordered_transp(l,j)             &
                                           * ao_coef_normalized_ordered_transp(m,i) 
    enddo
  enddo
 enddo
end


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



double precision function phi_j_erf_mu_r_phi(i,j,mu_in, C_center)
 implicit none
 BEGIN_DOC
! phi_j_erf_mu_r_phi  = int dr phi_j(r) [erf(mu  |r - C|)/|r-C|]  phi_i(r)
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: mu_in, C_center(3)
 integer :: num_A,power_A(3), num_b, power_B(3)
 double precision :: alpha, beta, A_center(3), B_center(3),contrib,NAI_pol_mult_erf
 integer :: n_pt_in,l,m
 n_pt_in = n_pt_max_integrals
 ! j 
 num_A = ao_nucl(j)
 power_A(1:3)= ao_power(j,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 ! i 
 num_B = ao_nucl(i)
 power_B(1:3)= ao_power(i,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)

 phi_j_erf_mu_r_phi = 0.d0
 do l=1,ao_prim_num(j)
  alpha = ao_expo_ordered_transp(l,j)
  do m=1,ao_prim_num(i)
    beta = ao_expo_ordered_transp(m,i)
    contrib = NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)  
    phi_j_erf_mu_r_phi += contrib *  ao_coef_normalized_ordered_transp(l,j)             &
                                  *  ao_coef_normalized_ordered_transp(m,i) 
  enddo
 enddo
end


subroutine erf_mu_gauss_xyz_ij_ao(i,j,mu, C_center, delta,gauss_ints)
 implicit none
 BEGIN_DOC
  ! gauss_ints(m) =   \int dr exp(-delta (r - C)^2 ) x/y/z * erf(mu |r-r'|)/ |r-r'| * AO_i(r') * AO_j(r')
  !
  ! with m = 1 ==> x, m = 2, m = 3 ==> z
  !
  !      m = 4 ==> no x/y/z
 END_DOC
 integer, intent(in) :: i,j
 double precision, intent(in) :: mu, C_center(3),delta
 double precision, intent(out):: gauss_ints(4)

 integer :: num_A,power_A(3), num_b, power_B(3)
 double precision :: alpha, beta, A_center(3), B_center(3),contrib,NAI_pol_mult_erf
 double precision :: xyz_ints(4)
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

 gauss_ints = 0.d0
 do l=1,ao_prim_num(j)
  alpha = ao_expo_ordered_transp(l,j)
  do m=1,ao_prim_num(i)
    beta = ao_expo_ordered_transp(m,i)
    call erf_mu_gauss_xyz(C_center,delta,mu,A_center,B_center,power_A,power_B,alpha,beta,n_pt_in,xyz_ints)
    do mm = 1, 4
     gauss_ints(mm) += xyz_ints(mm)  * ao_coef_normalized_ordered_transp(l,j)             &
                                     * ao_coef_normalized_ordered_transp(m,i) 
    enddo
  enddo
 enddo
end




subroutine erf_mu_gauss_xyz(D_center,delta,mu,A_center,B_center,power_A,power_B,alpha,beta,n_pt_in,xyz_ints)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) x * erf(mu |r-r'|)/ |r-r'| * (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  !   xyz_ints(1) = x , xyz_ints(2) = y, xyz_ints(3) = z, xyz_ints(4) = x^0 
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta,mu  ! pure gaussian "D" and mu parameter
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3),n_pt_in
 double precision, intent(out)   :: xyz_ints(4)

 double precision  :: NAI_pol_mult_erf
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr,contrib
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,mm
 integer           :: power_B_tmp(3)
 dim1=100
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 accu = 0.d0
 do lx = 0, iorder_a_new(1)
  coefx = A_new(lx,1)
  if(dabs(coefx).lt.thr)cycle
  iorder_tmp(1) = lx
  do ly = 0, iorder_a_new(2)
   coefy = A_new(ly,2)
   coefxy = coefx * coefy 
   if(dabs(coefxy).lt.thr)cycle
   iorder_tmp(2) = ly
   do lz = 0, iorder_a_new(3)
    coefz = A_new(lz,3)
    coefxyz = coefxy * coefz 
    if(dabs(coefxyz).lt.thr)cycle
    iorder_tmp(3) = lz
     power_B_tmp = power_B
     contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
     xyz_ints(4) += contrib * coefxyz ! usual term with no x/y/z 
                                      
     do mm = 1, 3 
      ! (x phi_i ) * phi_j 
      ! x * (x - B_x)^b_x = b_x (x - B_x)^b_x + 1 * (x - B_x)^{b_x+1}
      
      !
      ! first contribution :: B_x (x - B_x)^b_x :: usual integral multiplied by B_x
      power_B_tmp = power_B
      contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
      xyz_ints(mm) += contrib * B_center(mm) * coefxyz 
                                                       
      !
      ! second contribution :: (x - B_x)^(b_x+1) :: integral with b_x=>b_x+1 
      power_B_tmp(mm) += 1
      contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
      xyz_ints(mm) += contrib * coefxyz     
     enddo
   enddo
  enddo
 enddo
 xyz_ints *= fact_a_new 
end

