subroutine test_int
 implicit none
 BEGIN_DOC
! test the function NAI_pol_mult_erf_gauss_r12
 END_DOC
 include 'utils/constants.include.F'
 double precision :: A_center(3), B_center(3), C_center(3), D_center(3), r(3)
 double precision :: alpha,beta,delta,mu,NAI_pol_mult_erf_gauss_r12,overlap_gauss_r12
 double precision :: numerical,  analytical, weight,primitive_value_explicit
 double precision :: numerical_j,analytical_j,r_ij_bis
 double precision :: gauss_a, gauss_b, gauss_d, coulomb, r_ij,jastrow
 integer          :: power_A(3), power_B(3), power_D(3)
 integer          :: ipoint,ao_i,ao_j,num_A,num_B,i,j,k,l
 double precision :: abs_error_av, relat_error_av,icount
 double precision :: abs_error_j_av, relat_error_j_av
 ! C    :: center of the Coulomb 
 mu = 1.d0
 C_center = 0.d0
 C_center(2) =  0.1d0
 C_center(1) = -0.3d0
 ! D    :: center of gaussian "D"
!  D_center = C_center
 D_center = 0.d0
 D_center(1) = -0.4534d0
 D_center(3) =  0.8934d0
 ! delta :: exponent of gaussian "D" 
 delta = 0.5964d0
 ! power_D          :: == 0 exponent of the polynom for the gaussian "D" 
 power_D = 0

 abs_error_av   = 0.d0
 relat_error_av = 0.d0
 abs_error_j_av   = 0.d0
 relat_error_j_av = 0.d0
 icount = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
!  i = 4
!  j = 1 
   print*,'i,j',i,j
    num_A = ao_nucl(i)
    power_A(1:3)= ao_power(i,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
   
    num_B = ao_nucl(j)
    power_B(1:3)= ao_power(j,1:3)
    B_center(1:3) = nucl_coord(num_B,1:3)
    do l=1,ao_prim_num(i)
!     l = 1
     alpha = ao_expo_ordered_transp(l,i)     
     do k=1,ao_prim_num(i)
!      k = 1
      beta = ao_expo_ordered_transp(k,j)     
   
       ! analytical integral 
       analytical = NAI_pol_mult_erf_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
       analytical_j = overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
       ! numerical  integral 
       numerical   = 0.d0
       numerical_j = 0.d0
       icount += 1.d0
       do ipoint = 1, n_points_final_grid
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)
        r_ij = dsqrt( (C_center(1) - r(1))**2 + (C_center(2) - r(2))**2 + (C_center(3) - r(3))**2 )
        r_ij_bis = dsqrt( (D_center(1) - r(1))**2 + (D_center(2) - r(2))**2 + (D_center(3) - r(3))**2 )
        weight = final_weight_at_r_vector(ipoint)
        if(dabs(r_ij).lt.1.d-6)then
         coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_ij**2 / (3.d0 *sqpi) 
        else
         coulomb = derf(mu * r_ij)/r_ij
        endif
        jastrow = dexp(-delta * r_ij_bis * r_ij_bis)
        gauss_a = primitive_value_explicit(power_A,A_center,alpha,r)
        gauss_b = primitive_value_explicit(power_B,B_center,beta ,r)
        gauss_d = primitive_value_explicit(power_D,D_center,delta,r)
        numerical += weight * gauss_d * gauss_a * gauss_b * coulomb
        numerical_j+= weight * gauss_a * gauss_b * jastrow
       enddo
       abs_error_av += dabs(analytical - numerical)
       abs_error_j_av += dabs(analytical_j - numerical_j)
!!!!    print*,'numerical  = ',numerical
!!!!    print*,'analytical = ',analytical
       if(dabs(numerical).gt.1.d-10)then
        relat_error_av += dabs(analytical - numerical)/dabs(numerical)
        if(dabs(analytical - numerical)/dabs(numerical) .gt. 1.d-6 )then
         print*,'i,j',i,j
         print*,'l,k',l,k
         print*,'power_A = ',power_A
         print*,'power_B = ',power_B
         print*,'alpha, beta', alpha, beta
         print*,'numerical, analytical ',numerical,analytical
         print*,'error      = ',dabs(analytical - numerical),dabs(analytical - numerical)/dabs(numerical)
        endif
       endif 

       if(dabs(numerical_j).gt.1.d-10)then
        relat_error_j_av += dabs(analytical_j - numerical_j)/dabs(numerical_j)
        if(dabs(analytical_j - numerical_j)/dabs(numerical_j) .gt. 1.d-6 )then
         print*,'i,j',i,j
         print*,'l,k',l,k
         print*,'power_A = ',power_A
         print*,'power_B = ',power_B
         print*,'alpha, beta', alpha, beta
         print*,'numerical_j, analytical_j ',numerical_j,analytical_j
         print*,'error      = ',dabs(analytical_j - numerical_j),dabs(analytical_j - numerical_j)/dabs(numerical_j)
         stop
        endif
       endif 
     enddo  ! k
    enddo ! l
  enddo ! j 
 enddo ! i
 abs_error_av   = abs_error_av/icount
 relat_error_av = relat_error_av /icount
 print*,'abs_error_av     = ',abs_error_av
 print*,'relat_error_av   = ',relat_error_av
 abs_error_j_av   = abs_error_j_av/icount
 relat_error_j_av = relat_error_j_av /icount
 print*,'abs_error_j_av   = ',abs_error_j_av
 print*,'relat_error_j_av = ',relat_error_j_av

end

subroutine test_pol
 BEGIN_DOC
! test the function give_pol_in_r
 END_DOC
 implicit none
 include 'utils/constants.include.F'
 double precision :: A_center(3), B_center(3), P_center(3), r(3)
 double precision :: P_new(0:n_pt_max_integrals,3)
 double precision :: alpha,beta,alpha_new,fact_k
 double precision :: numerical, analytical, primitive_value_explicit, give_pol_in_r
 double precision :: gauss_a, gauss_b 
 integer          :: power_A(3), power_B(3), iorder(3)
 integer          :: ipoint
 ! A, B :: center of the gaussians "A", "B" 
 A_center = 0.d0
 B_center = 0.d0
 ! alpha/beta :: exponents of gaussians "A", "B" 
 !               ao  prim
 alpha = ao_expo(1 , 2   )
 beta  = ao_expo(2 , 1   )
 ! power_A, power_B :: exponents of the polynoms for gaussians "A", "B" 
 power_A = 0
 power_B = 0

 call give_explicit_poly_and_gaussian(P_new,P_center,alpha_new,fact_k,iorder,alpha,beta,power_A,power_B,A_center,B_center,n_pt_max_integrals)

 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  analytical = fact_k * give_pol_in_r(r,P_new,P_center, alpha_new,iorder, n_pt_max_integrals)
  gauss_a = primitive_value_explicit(power_A,A_center,alpha,r)
  gauss_b = primitive_value_explicit(power_B,B_center,beta ,r)
  numerical = gauss_b * gauss_a
  if(dabs(numerical).gt.1.d-10)then
   if(dabs(numerical - analytical).gt.1.d-10)then
    print*,numerical,analytical,fact_k
   endif
  endif
 enddo


end

subroutine test_int_bis
 implicit none
 BEGIN_DOC
! test the function erf_mu_stg_gauss_int_phi_ij ovlp_stg_gauss_int_phi_ij
 END_DOC
 include 'utils/constants.include.F'
 double precision :: A_center(3), B_center(3), alpha,beta      ! basis functions
 integer          :: power_A(3), power_B(3)                    ! basis functions 
 double precision :: D_center(3), gama,delta  ! Jastrow factor 
 double precision :: C_center(3), mu          ! Coulomb erf(mu*x)/x

 double precision :: ovlp_stg_gauss_int_phi_ij,erf_mu_stg_gauss_int_phi_ij
 double precision :: numerical,  analytical, weight,primitive_value_explicit
 double precision :: numerical_c,analytical_c,r_ij_bis,r(3)
 double precision :: gauss_a, gauss_b, gauss_d, coulomb, r_ij,jastrow
 integer          :: ipoint,ao_i,ao_j,num_A,num_B,i,j,k,l
 double precision :: abs_error_av, relat_error_av,icount
 double precision :: abs_error_c_av, relat_error_c_av,slater_fit_gam
 ! C    :: center of the Coulomb 
 mu = 1.d0
 C_center = 0.d0
 C_center(:) = nucl_coord_transp(:,1)
! C_center(2) =  0.1d0
! C_center(1) = -0.3d0
 ! D    :: center of jastrow factor 
  D_center = C_center
! D_center = 0.d0
! D_center(1) = -0.2d0
! D_center(3) =  0.2d0 
 gama = 1.2d0 ! exponent of the Slater of the Jastrow factor 
 delta = 1.0 ! exponent of the gaussian in the Jastrow factor 

 abs_error_av   = 0.d0
 relat_error_av = 0.d0
 abs_error_c_av   = 0.d0
 relat_error_c_av = 0.d0
 icount = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   print*,'i,j',i,j
    num_A = ao_nucl(i)
    power_A(1:3)= ao_power(i,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
   
    num_B = ao_nucl(j)
    power_B(1:3)= ao_power(j,1:3)
    B_center(1:3) = nucl_coord(num_B,1:3)
    do l=1,ao_prim_num(i)
!     l = 1
     alpha = ao_expo_ordered_transp(l,i)     
     do k=1,ao_prim_num(j)
!      k = 1
      beta = ao_expo_ordered_transp(k,j)     
   
       ! analytical integral 
       analytical = ovlp_stg_gauss_int_phi_ij(D_center,gama,delta,A_center,B_center,power_A,power_B,alpha,beta)
       analytical_c = erf_mu_stg_gauss_int_phi_ij(D_center,gama,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)

       ! numerical  integral 
       numerical   = 0.d0
       numerical_c = 0.d0
       icount += 1.d0
       do ipoint = 1, n_points_final_grid
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)
        weight = final_weight_at_r_vector(ipoint)
        ! distance from r to D
        r_ij     = dsqrt( (D_center(1) - r(1))**2 + (D_center(2) - r(2))**2 + (D_center(3) - r(3))**2 )
        ! distance from r to C
        r_ij_bis = dsqrt( (C_center(1) - r(1))**2 + (C_center(2) - r(2))**2 + (C_center(3) - r(3))**2 )
        if(dabs(r_ij_bis).lt.1.d-6)then
         coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_ij_bis**2 / (3.d0 *sqpi) 
        else
         coulomb = derf(mu * r_ij_bis)/r_ij_bis
        endif
        ! jastrow factor 
        jastrow = slater_fit_gam(r_ij,gama) * dexp(-delta * r_ij * r_ij)
        gauss_a = primitive_value_explicit(power_A,A_center,alpha,r)
        gauss_b = primitive_value_explicit(power_B,B_center,beta ,r)
        numerical += weight * jastrow * gauss_a * gauss_b 
        numerical_c += weight * jastrow * gauss_a * gauss_b * coulomb
       enddo
       abs_error_av += dabs(analytical - numerical)
       abs_error_c_av += dabs(analytical_c - numerical_c)
!       abs_error_j_av += dabs(analytical_j - numerical_j)
       if(dabs(numerical).gt.1.d-10)then
        relat_error_av += dabs(analytical - numerical)/dabs(numerical)
        if(dabs(analytical - numerical)/dabs(numerical) .gt. 1.d-6 )then
         print*,'i,j',i,j
         print*,'l,k',l,k
         print*,'power_A = ',power_A
         print*,'power_B = ',power_B
         print*,'alpha, beta', alpha, beta
         print*,'numerical, analytical ',numerical,analytical
         print*,'error      = ',dabs(analytical - numerical),dabs(analytical - numerical)/dabs(numerical)
        endif
       endif 

       if(dabs(numerical_c).gt.1.d-10)then
        relat_error_c_av += dabs(analytical_c - numerical_c)/dabs(numerical_c)
        if(dabs(analytical_c - numerical_c)/dabs(numerical_c) .gt. 1.d-6 )then
         print*,'i,j',i,j
         print*,'l,k',l,k
         print*,'power_A = ',power_A
         print*,'power_B = ',power_B
         print*,'alpha, beta', alpha, beta
         print*,'numerical_c, analytical_c ',numerical_c,analytical_c
         print*,'error      = ',dabs(analytical_c - numerical_c),dabs(analytical_c - numerical_c)/dabs(numerical_c)
!         stop
        endif
       endif 

     enddo  ! k
    enddo ! l
  enddo ! j 
 enddo ! i
 abs_error_av   = abs_error_av/icount
 relat_error_av = relat_error_av /icount
 print*,'abs_error_av     = ',abs_error_av
 print*,'relat_error_av   = ',relat_error_av
 abs_error_c_av   = abs_error_c_av/icount
 relat_error_c_av = relat_error_c_av /icount
 print*,'abs_error_c_av   = ',abs_error_c_av
 print*,'relat_error_c_av = ',relat_error_c_av

end

subroutine test_int_full_j
 implicit none
 BEGIN_DOC
! test the function erf_mu_stg_gauss_int_phi_ij ovlp_stg_gauss_int_phi_ij
 END_DOC
 include 'utils/constants.include.F'
 double precision :: A_center(3), B_center(3), alpha,beta      ! basis functions
 integer          :: power_A(3), power_B(3)                    ! basis functions 
 double precision :: D_center(3), gama,delta  ! Jastrow factor 
 double precision :: C_center(3), mu          ! Coulomb erf(mu*x)/x

 double precision :: ovlp_stg_gauss_int_phi_ij,erf_mu_stg_gauss_int_phi_ij
 double precision :: numerical,  analytical, weight,primitive_value_explicit
 double precision :: numerical_c,analytical_c,r_ij_bis,r(3)
 double precision :: gauss_a, gauss_b, gauss_d, coulomb, r_ij,jastrow
 integer          :: ipoint,ao_i,ao_j,num_A,num_B,i,j,k,l
 double precision :: abs_error_av, relat_error_av,icount
 double precision :: abs_error_c_av, relat_error_c_av,slater_fit_gam
 double precision :: full_jastrow,exp_dl
 integer :: n_taylor
 n_taylor = 10
 double precision, allocatable :: array_ints_overlap(:),array_ints_erf(:)
 allocate(array_ints_overlap(0:n_taylor),array_ints_erf(0:n_taylor))
 ! C    :: center of the Coulomb 
 mu = 1.d0
! C_center(:) = nucl_coord_transp(:,1)
 C_center = 0.d0
 C_center(2) =  0.1d0
 C_center(1) = -0.3d0
 ! D    :: center of jastrow factor 
!  D_center = C_center
 D_center = 0.d0
 D_center(1) = -0.2d0
 D_center(3) =  0.2d0 
 gama = 1.2d0 ! exponent of the Slater of the Jastrow factor 
 delta = 1.0 ! exponent of the gaussian in the Jastrow factor 

 abs_error_av   = 0.d0
 relat_error_av = 0.d0
 abs_error_c_av   = 0.d0
 relat_error_c_av = 0.d0
 icount = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   print*,'i,j',i,j
    num_A = ao_nucl(i)
    power_A(1:3)= ao_power(i,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)
   
    num_B = ao_nucl(j)
    power_B(1:3)= ao_power(j,1:3)
    B_center(1:3) = nucl_coord(num_B,1:3)
    do l=1,ao_prim_num(i)
!     l = 1
     alpha = ao_expo_ordered_transp(l,i)     
     do k=1,ao_prim_num(j)
!      k = 1
      beta = ao_expo_ordered_transp(k,j)     
   
       ! analytical integral 
       call exp_dl_ovlp_stg_phi_ij(D_center,gama,delta,A_center,B_center,power_A,power_B,alpha,beta,n_taylor,array_ints_overlap,analytical)
       call exp_dl_erf_stg_phi_ij( D_center,gama,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu,n_taylor,array_ints_erf,analytical_c)

       ! numerical  integral 
       numerical   = 0.d0
       numerical_c = 0.d0
       icount += 1.d0
       do ipoint = 1, n_points_final_grid
        r(1) = final_grid_points(1,ipoint)
        r(2) = final_grid_points(2,ipoint)
        r(3) = final_grid_points(3,ipoint)
        weight = final_weight_at_r_vector(ipoint)
        ! distance from r to D
        r_ij     = dsqrt( (D_center(1) - r(1))**2 + (D_center(2) - r(2))**2 + (D_center(3) - r(3))**2 )
        ! distance from r to C
        r_ij_bis = dsqrt( (C_center(1) - r(1))**2 + (C_center(2) - r(2))**2 + (C_center(3) - r(3))**2 )
        if(dabs(r_ij_bis).lt.1.d-6)then
         coulomb = 2.d0 * mu / sqpi - 2.d0 * mu**3 * r_ij_bis**2 / (3.d0 *sqpi) 
        else
         coulomb = derf(mu * r_ij_bis)/r_ij_bis
        endif
        ! jastrow factor 
        jastrow = slater_fit_gam(r_ij,gama) * dexp(-delta * r_ij * r_ij)
        full_jastrow = exp_dl(-jastrow, n_taylor)
!        full_jastrow = dexp(-jastrow)
        gauss_a = primitive_value_explicit(power_A,A_center,alpha,r)
        gauss_b = primitive_value_explicit(power_B,B_center,beta ,r)
        numerical += weight * full_jastrow * gauss_a * gauss_b 
        numerical_c += weight * full_jastrow * gauss_a * gauss_b * coulomb
       enddo
       abs_error_av += dabs(analytical - numerical)
       if(dabs(numerical).gt.1.d-10)then
        relat_error_av += dabs(analytical - numerical)/dabs(numerical)
        if(dabs(analytical - numerical)/dabs(numerical) .gt. 1.d-6 )then
         print*,'i,j',i,j
         print*,'l,k',l,k
         print*,'power_A = ',power_A
         print*,'power_B = ',power_B
         print*,'alpha, beta', alpha, beta
         print*,'numerical, analytical ',numerical,analytical
         print*,'error      = ',dabs(analytical - numerical),dabs(analytical - numerical)/dabs(numerical)
        endif
       endif 

       abs_error_c_av += dabs(analytical_c - numerical_c)
       if(dabs(numerical_c).gt.1.d-10)then
        relat_error_c_av += dabs(analytical_c - numerical_c)/dabs(numerical_c)
        if(dabs(analytical_c - numerical_c)/dabs(numerical_c) .gt. 1.d-6 )then
         print*,'i,j',i,j
         print*,'l,k',l,k
         print*,'power_A = ',power_A
         print*,'power_B = ',power_B
         print*,'alpha, beta', alpha, beta
         print*,'numerical_c, analytical_c ',numerical_c,analytical_c
         print*,'error      = ',dabs(analytical_c - numerical_c),dabs(analytical_c - numerical_c)/dabs(numerical_c)
        endif
       endif 

     enddo  ! k
    enddo ! l
  enddo ! j 
 enddo ! i
 abs_error_av   = abs_error_av/icount
 relat_error_av = relat_error_av /icount
 print*,'abs_error_av     = ',abs_error_av
 print*,'relat_error_av   = ',relat_error_av
 abs_error_c_av   = abs_error_c_av/icount
 relat_error_c_av = relat_error_c_av /icount
 print*,'abs_error_c_av   = ',abs_error_c_av
 print*,'relat_error_c_av = ',relat_error_c_av

end
