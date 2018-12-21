!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function correction_to_on_top_from_UEG(mu,r,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: mu,r(3)
 double precision :: rho_a(N_states),rho_b(N_states)
 double precision :: g0_UEG_mu_inf, g0_UEG_mu 
 call dm_dft_alpha_beta_at_r(r,rho_a,rho_b)

 correction_to_on_top_from_UEG = g0_UEG_mu_inf(rho_a(istate),rho_b(istate)) / g0_UEG_mu(mu,rho_a(istate),rho_b(istate)) 

end




double precision function g0_UEG_mu_inf(rho_a,rho_b)
 implicit none
 double precision, intent(in) :: rho_a,rho_b
 double precision :: rho,pi,x
 double precision :: B, C, D, E, d2, rs, ahd
 rho = rho_a+rho_b
 pi = 4d0 * datan(1d0)
 ahd = -0.36583d0
 d2 = 0.7524d0
 B = -2d0 * ahd - d2
 C = 0.08193d0
 D = -0.01277d0
 E = 0.001859d0
 rs = (3d0 / 4d0*pi*rho)**(1d0/3d0)
 x = -d2*rs

 g0_UEG_mu_inf= 0.5d0 * (1d0- B*rs + C*rs**2d0 + D*rs**3d0 + E*rs**4d0)*exp(x)
end



double precision function g0_UEG_mu(mu,rho_a,rho_b)
 implicit none
 double precision, intent(in) :: rho_a,rho_b,mu
 double precision :: zeta,pi,rho,x,alpha
 double precision :: B, C, D, E, d2, rs, ahd, h, kf
 pi = 4d0 * datan(1d0)
 rho = rho_a+rho_b
!zeta = (rho_a-rho_b)/(rho_a+rho_b)
 alpha = (4d0/(9d0*pi))**(1d0/3d0)
 ahd = -0.36583d0
 d2 = 0.7524d0
 B = -2d0 * ahd - d2
 C = 0.08193d0
 D = -0.01277d0
 E = 0.001859d0
 rs = (3d0 / 4d0*pi*rho)**(1d0/3d0)
 kf = (alpha*rs)**(-1d0)
 zeta = mu / kf
 x = -d2*rs*h(zeta)/ahd 
 g0_UEG_mu = (exp(x)/2d0) * (1d0- B*(h(zeta)/ahd)*rs + C*((h(zeta)**2d0)/(ahd**2d0))*(rs**2d0) + D*((h(zeta)**3d0)/(ahd**3d0))*(rs**3d0) + E*((h(zeta)**4d0)/(ahd**4d0))*(rs**4d0) )
 
end



double precision function h(zeta)
 implicit none
 double precision, intent(in) :: zeta
 double precision :: pi
 double precision :: a1, a2, b1, b2, b3, ahd, alpha
 pi = 4d0 * datan(1d0)
 ahd = -0.36583d0
 alpha = (4d0/(9d0*pi))**(1d0/3d0)
 a1 = -(6d0*alpha/pi)*(1d0-log(2d0))
 b1 = 1.4919d0
 b3 = 1.91528d0
 a2 = ahd * b3
 b2 = (a1 - (b3*alpha/sqrt(pi)))/ahd

 h = (a1*zeta**2d0 + a2*zeta**3d0) / (1d0 + b1*zeta + b2*zeta**2d0 + b3*zeta**3d0)
end



