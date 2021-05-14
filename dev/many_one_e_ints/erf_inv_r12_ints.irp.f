
double precision function NAI_pol_mult_erf_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !   \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: C_center(3),mu      ! coulomb center "C" and "mu" in the erf(mu*x)/x function
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 double precision  :: NAI_pol_mult_erf
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3)
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
    accu += coefxyz * NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B,alpha_new,beta,C_center,n_pt_max_integrals,mu)
   enddo
  enddo
 enddo
 NAI_pol_mult_erf_gauss_r12 = fact_a_new * accu 
end


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



subroutine erfc_mu_gauss_xyz_ij_ao(i,j,mu, C_center, delta,gauss_ints)
 implicit none
 BEGIN_DOC
  ! gauss_ints(m) =   \int dr exp(-delta (r - C)^2 ) x/y/z * ( 1 - erf(mu |r-r'|))/ |r-r'| * AO_i(r') * AO_j(r')
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
    call erfc_mu_gauss_xyz(C_center,delta,mu,A_center,B_center,power_A,power_B,alpha,beta,n_pt_in,xyz_ints)
    do mm = 1, 4
     gauss_ints(mm) += xyz_ints(mm)  * ao_coef_normalized_ordered_transp(l,j)             &
                                     * ao_coef_normalized_ordered_transp(m,i) 
    enddo
  enddo
 enddo
end




subroutine erfc_mu_gauss_xyz(D_center,delta,mu,A_center,B_center,power_A,power_B,alpha,beta,n_pt_in,xyz_ints)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) x * (1 - erf(mu |r-r'|))/ |r-r'| * (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
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
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr,contrib,contrib_inf,mu_inf
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,mm
 integer           :: power_B_tmp(3)
 dim1=100
 mu_inf = 1.d+10
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 xyz_ints = 0.d0
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
     contrib_inf = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu_inf)  
     xyz_ints(4) += (contrib_inf - contrib) * coefxyz ! usual term with no x/y/z 
                                      
     do mm = 1, 3 
      ! (x phi_i ) * phi_j 
      ! x * (x - B_x)^b_x = B_x (x - B_x)^b_x + 1 * (x - B_x)^{b_x+1}
      
      !
      ! first contribution :: B_x (x - B_x)^b_x :: usual integral multiplied by B_x
      power_B_tmp = power_B
      contrib_inf = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu_inf)  
      contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
      xyz_ints(mm) += (contrib_inf - contrib) * B_center(mm) * coefxyz 
                                                       
      !
      ! second contribution :: (x - B_x)^(b_x+1) :: integral with b_x=>b_x+1 
      power_B_tmp(mm) += 1
      contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
      contrib_inf = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu_inf)  
      xyz_ints(mm) += (contrib_inf -  contrib) * coefxyz     
     enddo
   enddo
  enddo
 enddo
 xyz_ints *= fact_a_new 
end

subroutine NAI_pol_x_mult_erf_ao(i_ao,j_ao,mu_in,C_center,ints)
 implicit none
  BEGIN_DOC
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr y * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  !
  ! $\int_{-\infty}^{infty} dr z * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  END_DOC
 include 'utils/constants.include.F'                                                                                                                                  
 integer, intent(in) :: i_ao,j_ao
 double precision, intent(in) :: mu_in, C_center(3)
 double precision, intent(out):: ints(3)
 double precision               :: A_center(3), B_center(3),integral, alpha,beta
 double precision               :: NAI_pol_mult_erf
 integer                        :: i,j,num_A,num_B, power_A(3), power_B(3), n_pt_in, power_xA(3),m
 num_A = ao_nucl(i_ao)
 power_A(1:3)= ao_power(i_ao,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j_ao)
 power_B(1:3)= ao_power(j_ao,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 n_pt_in = n_pt_max_integrals

 ints = 0.d0

 do i = 1, ao_prim_num(i_ao)
  alpha = ao_expo_ordered_transp(i,i_ao)
   do m = 1, 3
    power_xA = power_A
    ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
    power_xA(m) += 1
    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      ! First term = (x-Ax)**(ax+1)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_xA,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints(m) += integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
      ! Second term = Ax * (x-Ax)**(ax)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints(m) += A_center(m) * integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
    enddo
  enddo
 enddo
end

subroutine NAI_pol_x_specify_mult_erf_ao(i_ao,j_ao,mu_in,C_center,m,ints)
 implicit none
  BEGIN_DOC
  ! Computes the following integral :
  ! $\int_{-\infty}^{infty} dr x * \chi_i(r) \chi_j(r) \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$.
  END_DOC
 include 'utils/constants.include.F'                                                                                                                                  
 integer, intent(in) :: i_ao,j_ao,m
 double precision, intent(in) :: mu_in, C_center(3)
 double precision, intent(out):: ints
 double precision               :: A_center(3), B_center(3),integral, alpha,beta
 double precision               :: NAI_pol_mult_erf
 integer                        :: i,j,num_A,num_B, power_A(3), power_B(3), n_pt_in, power_xA(3)
 num_A = ao_nucl(i_ao)
 power_A(1:3)= ao_power(i_ao,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j_ao)
 power_B(1:3)= ao_power(j_ao,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 n_pt_in = n_pt_max_integrals

 ints = 0.d0

 do i = 1, ao_prim_num(i_ao)
  alpha = ao_expo_ordered_transp(i,i_ao)
    power_xA = power_A
    ! x * phi_i(r) = x * (x-Ax)**ax = (x-Ax)**(ax+1) + Ax * (x-Ax)**ax
    power_xA(m) += 1
    do j = 1, ao_prim_num(j_ao)
      beta = ao_expo_ordered_transp(j,j_ao)
      ! First term = (x-Ax)**(ax+1)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_xA,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints += integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
      ! Second term = Ax * (x-Ax)**(ax)
      integral =  NAI_pol_mult_erf(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in,mu_in)
      ints += A_center(m) * integral * ao_coef_normalized_ordered_transp(j,j_ao)*ao_coef_normalized_ordered_transp(i,i_ao)
    enddo
 enddo
end

