
 BEGIN_PROVIDER [double precision, mu_of_r_basis_hf, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, grad_mu_of_r_basis_hf, (3,n_points_final_grid) ]
 implicit none 
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals 
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,r(3),mu,mu_tmp,mu_basis_hf,damped_mu,mu_min,grad_mu(3),dx
 PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals damped_mu_of_r mu_of_r_min
 print*,'providing mu_of_r_basis_hf ...'
 call wall_time(wall0)
 mu_min = mu_erf
 dx = 1.d-5
! !$OMP PARALLEL DO &
! !$OMP DEFAULT (NONE)  &
! !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf,r,mu,mu_tmp,grad_mu,dx) & 
! !$OMP SHARED (n_points_final_grid,mu_of_r_basis_hf,final_grid_points,mu_of_r_min,damped_mu_of_r,mu_min, grad_mu_of_r_basis_hf) 
 do ipoint = 1, n_points_final_grid
  r(:) = final_grid_points(:,ipoint)
  mu_tmp = mu_basis_hf(r)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
  else
   mu = max(mu_tmp,mu_of_r_min)
  endif
  mu_of_r_basis_hf(ipoint) = mu 

  if(damped_mu_of_r)then
   call get_grad_damped_mu_hf(r,dx,mu_min,grad_mu)
  else
   call get_grad_mu_hf(r,dx,grad_mu)
  endif
  grad_mu_of_r_basis_hf(:,ipoint) = grad_mu(:)
 enddo
! !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_basis_hf = ',wall1-wall0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_of_r_extra_grid_basis_hf, (n_points_extra_final_grid) ]
 implicit none 
 BEGIN_DOC
 ! mu(r) computed with a HF wave function (assumes that HF MOs are stored in the EZFIO)
 !
 ! corresponds to Eq. (37) of J. Chem. Phys. 149, 194301 (2018) but for \Psi^B = HF^B
 !
 ! !!!!!! WARNING !!!!!! if no_core_density == .True. then all contributions from the core orbitals 
 !
 ! in the two-body density matrix are excluded
 END_DOC
 integer :: ipoint
 double precision :: wall0,wall1,f_hf,on_top,w_hf,r(3),mu,mu_tmp,mu_basis_hf,damped_mu,mu_min
 PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals damped_mu_of_r mu_of_r_min
 print*,'providing mu_of_r_extra_grid_basis_hf ...'
 call wall_time(wall0)
 mu_min = mu_erf
! !$OMP PARALLEL DO &
! !$OMP DEFAULT (NONE)  &
! !$OMP PRIVATE (ipoint,f_hf,on_top,w_hf,r,mu,mu_tmp) & 
! !$OMP ShARED (n_points_extra_final_grid,mu_of_r_extra_grid_basis_hf,final_grid_points_extra,mu_of_r_min,damped_mu_of_r,mu_min) 
 do ipoint = 1, n_points_extra_final_grid
  r(:) = final_grid_points_extra(:,ipoint)
  mu_tmp = mu_basis_hf(r)
  if(damped_mu_of_r)then
   mu = damped_mu(mu_tmp,mu_min)
  else
   mu = max(mu_tmp,mu_of_r_min)
  endif
  mu_of_r_extra_grid_basis_hf(ipoint) = mu 
 enddo
! !$OMP END PARALLEL DO
 call wall_time(wall1)
 print*,'Time to provide mu_of_r_extra_grid_basis_hf = ',wall1-wall0
 END_PROVIDER 

double precision function mu_basis_damped(r,mu_min)
 implicit none
 double precision, intent(in) :: r(3), mu_min
 double precision :: mu_basis_hf, damped_mu, mu_tmp
 mu_tmp = mu_basis_hf(r)
 mu_basis_damped = damped_mu(mu_tmp,mu_min)
end

double precision function mu_basis_hf(r)
 implicit none
 double precision, intent(in) :: r(3)
 double precision :: f_hf, on_top, w_hf
 include 'constants.include.F'
 call f_HF_valence_ab(r,r,f_hf,on_top)
 if(on_top.le.1.d-12.or.f_hf.le.0.d0.or.f_hf * on_top.lt.0.d0)then
   w_hf   = 1.d+10
 else 
   w_hf  = f_hf /  on_top
 endif
 mu_basis_hf =  w_hf * sqpi * 0.5d0
end

subroutine get_grad_mu_hf(r,dx,grad_mu)
  implicit none
  double precision, intent(in) :: r(3),dx
  double precision, intent(out):: grad_mu(3)
  double precision :: r1(3),mu_plus,mu_minus,mu_basis_hf
  integer :: m
  do m = 1, 3 ! compute grad mu
   r1 = r
   r1(m) += dx 
   mu_plus = mu_basis_hf(r1)
   r1 = r 
   r1(m) -= dx 
   mu_minus = mu_basis_hf(r1)
   grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo
end


subroutine get_grad_damped_mu_hf(r,dx,mu_min,grad_mu)
  implicit none
  double precision, intent(in) :: r(3),dx,mu_min
  double precision, intent(out):: grad_mu(3)
  double precision :: r1(3),mu_plus,mu_minus,mu_basis_damped
  integer :: m
  do m = 1, 3 ! compute grad mu
   r1 = r
   r1(m) += dx 
   mu_plus = mu_basis_damped(r1,mu_min)
   r1 = r 
   r1(m) -= dx 
   mu_minus = mu_basis_damped(r1,mu_min)
   grad_mu(m) = (mu_plus - mu_minus)/(2.d0 * dx)
  enddo
end

