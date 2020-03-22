 BEGIN_PROVIDER [double precision, Energy_c_md_on_top, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density.
  ! Ec_md_on_top = (alpha/mu**3) * int n2(r,r) dr  where alpha = sqrt(2pi)*(-2+sqrt(2)) 
 END_DOC
 implicit none 
 integer :: istate
 double precision :: pi,mu
 mu = mu_erf_dft
 pi = 4.d0 * datan(1.d0)
 Energy_c_md_on_top = ((-2.d0+sqrt(2.d0))*sqrt(2.d0*pi)/(3.d0*(mu**3)))*integral_on_top/(1.d0 + 2.d0 / (dsqrt(pi)*mu ))
 END_PROVIDER


 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_spin_pol]
 BEGIN_DOC
! new ECMD functional for spin polarized electrons
 END_DOC
 implicit none
 integer :: i_point
 double precision :: weight
 double precision :: constant
 constant = (-3.d0 + 2.d0 * dsqrt(2.d0)) * dsqrt(dacos(-1.d0)) / (10.d0*dsqrt(2.d0) * mu_erf_dft**5)
 Energy_c_md_on_top_spin_pol = 0.d0
 do i_point = 1, n_points_final_grid
  weight = final_weight_at_r_vector(i_point)
  Energy_c_md_on_top_spin_pol += nabla2_n2_hf_alpha_alpha(i_point) * weight 
 enddo
 ! extrapolation of the on top
 Energy_c_md_on_top_spin_pol *= constant / (1.d0 + 2.d0 / (3.d0 * dsqrt(dacos(-1.d0))*mu_erf_dft ))
 END_PROVIDER 

