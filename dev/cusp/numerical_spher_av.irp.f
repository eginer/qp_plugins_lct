double precision function spher_av_n2_hf(r1,r12)
 BEGIN_DOC
! computes the spherical average of n(r1)n(r2) for a sphere of radius r12
 END_DOC
 implicit none
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: r2(3),on_top_00
 spher_av_n2_hf = 0.d0

 do k = 1, n_points_integration_angular
  r2(1)= r1(1)+ r12* angular_quadrature_points(k,1)
  r2(2)= r1(2)+ r12* angular_quadrature_points(k,2)
  r2(3)= r1(3)+ r12* angular_quadrature_points(k,3)

  call two_body_dm_hf_r1r2(r1,r2,on_top_00)
  spher_av_n2_hf += on_top_00 * weights_angular_points(k)
 enddo
 spher_av_n2_hf = 0.25d0/dacos(-1.d0) * spher_av_n2_hf

end 

double precision function spher_av_n2_hartree(r1,r12)
 BEGIN_DOC
! computes the spherical average of n(r1)n(r2) for a sphere of radius r12
 END_DOC
 implicit none
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: r2(3),rhoa_r1,rhob_r1,rho_tot_r1
 double precision :: rhoa_r2,rhob_r2,rho_tot_r2
 spher_av_n2_hartree = 0.d0

 call dm_dft_alpha_beta_at_r(r1,rhoa_r1,rhob_r1)
 rho_tot_r1 = rhoa_r1 + rhob_r1
 do k = 1, n_points_integration_angular
  r2(1)= r1(1)+ r12* angular_quadrature_points(k,1)
  r2(2)= r1(2)+ r12* angular_quadrature_points(k,2)
  r2(3)= r1(3)+ r12* angular_quadrature_points(k,3)

  call dm_dft_alpha_beta_at_r(r2,rhoa_r2,rhob_r2)
  rho_tot_r2 = rhoa_r2 + rhob_r2

  spher_av_n2_hartree += rho_tot_r2 * rho_tot_r1 * weights_angular_points(k)
 enddo
 spher_av_n2_hartree = 0.25d0/dacos(-1.d0) * spher_av_n2_hartree

end 

subroutine spher_av_n2_psi(r1,r12,spher_av_n2,spher_av_j_r12)
 BEGIN_DOC
! computes the spherical average of n2(r1,r2) for a sphere of radius r12
 END_DOC
 implicit none
 double precision, intent(in)  :: r1(3),r12
 double precision, intent(out) :: spher_av_n2,spher_av_j_r12
 integer :: k
 double precision :: r2(3),on_top_00,on_top_hf
 spher_av_n2= 0.d0
 spher_av_j_r12 = 0.d0
 do k = 1, n_points_integration_angular
  r2(1)= r1(1)+ r12* angular_quadrature_points(k,1)
  r2(2)= r1(2)+ r12* angular_quadrature_points(k,2)
  r2(3)= r1(3)+ r12* angular_quadrature_points(k,3)

  call two_body_r1r2(r1,r2,1,on_top_00)
  call two_body_dm_hf_r1r2(r1,r2,on_top_hf)

  spher_av_n2    += on_top_00           * weights_angular_points(k)
  spher_av_j_r12 += on_top_00/on_top_hf * weights_angular_points(k)
 enddo
 spher_av_n2    = 0.25d0/dacos(-1.d0) * spher_av_n2
 spher_av_j_r12 = 0.25d0/dacos(-1.d0) * spher_av_j_r12

end 

