program plot_on_top
 implicit none
 io_two_body_rdm_ab = "Read"
 touch io_two_body_rdm_ab
 call routine
end
subroutine routine
 implicit none
 integer :: nx,i
 double precision :: dx,r(3),xmax,xmin
 double precision :: mu,f_psi,on_top, on_top_extrap, g0_UEG_mu_inf
 double precision :: mu_correction_of_on_top, rho_a, rho_b, on_top_ueg
 integer :: istate 
 integer                        :: i_unit_output,getUnitAndOpen
 character*(128)                :: output
 PROVIDE ezfio_filename
 output=trim(ezfio_filename)//'.on_top_psi'
 i_unit_output = getUnitAndOpen(output,'w')

 istate = 1
 nx = 1000
 xmax =  3.d0
 xmin = -xmax
 dx = (xmax - xmin)/dble(nx)
 
 r = 0.d0 
 r(3) = xmin 
 do i = 1, nx 
  call give_mu_of_r_cas(r,istate,mu,f_psi,on_top)
  ! on_top of the wave function 
  call dm_dft_alpha_beta_at_r(r,rho_a,rho_b)

! Extrapolated on-top pair density 
  on_top_extrap = 2.d0 * mu_correction_of_on_top(mu,on_top)
! We take the on-top pair density of the UEG which is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
  on_top_ueg = 4.d0 * rho_a * rho_b * g0_UEG_mu_inf(rho_a,rho_b)
  write(i_unit_output,'(X,100(F16.10,X))')r(3),2.d0 * on_top, on_top_ueg, on_top_extrap, rho_a, rho_b
  r(3) += dx
 enddo
end
