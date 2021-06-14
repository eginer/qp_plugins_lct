program plot_mu_of_r_tc
 implicit none
 double precision :: x,r(3),dx,xmin,xmax,r_tmp(3)
 double precision :: mu_lda,mu_minus,mu_plus,grad_mu,rho_a_hf,rho_b_hf,mu
 double precision :: mu_lda_damped,mu_min
 integer :: i,nx,m
 nx = 1000
 xmax = 10.d0
 xmin = -10.d0
 dx = (xmax - xmin)/dble(nx)
 x = xmin 
 mu_min = 0.5d0
 do i = 1, nx
  r = 0.d0
  r(3) = x
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  mu = mu_lda(rho_a_hf,rho_b_hf)
  grad_mu = 0.d0
  do m = 1, 3
   r_tmp = r
   r_tmp(m) += dx 
   call dm_dft_alpha_beta_at_r(r_tmp,rho_a_hf,rho_b_hf)
   mu_plus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
   r_tmp = r
   r_tmp(m) -= dx 
   call dm_dft_alpha_beta_at_r(r_tmp,rho_a_hf,rho_b_hf)
   mu_minus = mu_lda_damped(rho_a_hf,rho_b_hf,mu_min)
   grad_mu += ((mu_plus - mu_minus)/(2.d0 * dx))**2
  enddo
  call dm_dft_alpha_beta_at_r(r,rho_a_hf,rho_b_hf)
  write(33,'(100(F16.10,X))')x,rho_a_hf + rho_b_hf, mu, mu_lda_damped(rho_a_hf,rho_b_hf,mu_min),dsqrt(grad_mu)
  x += dx
 enddo

end

