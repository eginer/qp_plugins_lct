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

double precision function overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 double precision  :: overlap_x,overlap_y,overlap_z,overlap
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1
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
    call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
    accu += coefxyz * overlap
   enddo
  enddo
 enddo
 overlap_gauss_r12 = fact_a_new * accu 
end


double precision function overlap_gauss_r12_ao(D_center,delta,i,j)
 implicit none
 integer, intent(in) :: i,j
 double precision, intent(in) :: D_center(3), delta

 integer :: num_a,num_b,power_A(3), power_B(3),l,k
 double precision :: A_center(3), B_center(3),overlap_gauss_r12,alpha,beta,analytical_j
 num_A = ao_nucl(i)
 power_A(1:3)= ao_power(i,1:3)
 A_center(1:3) = nucl_coord(num_A,1:3)
 num_B = ao_nucl(j)
 power_B(1:3)= ao_power(j,1:3)
 B_center(1:3) = nucl_coord(num_B,1:3)
 overlap_gauss_r12_ao = 0.d0
 do l=1,ao_prim_num(i)
  alpha = ao_expo_ordered_transp(l,i)     
  do k=1,ao_prim_num(j)
   beta = ao_expo_ordered_transp(k,j)     
   analytical_j = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
   overlap_gauss_r12_ao += analytical_j *  ao_coef_normalized_ordered_transp(l,i)             &
                                        *  ao_coef_normalized_ordered_transp(k,j) 
  enddo 
 enddo
end
