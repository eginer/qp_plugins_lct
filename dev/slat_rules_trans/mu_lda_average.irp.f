program test_density
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
!  read_wf = .True.
!  touch read_wf
!  read_rl_eigv = .True.
!  touch read_rl_eigv
!  call routine_write
   print*,'average_mu_lda = ',average_mu_lda
   print*,'average_mu_lda * 1/2 = ',average_mu_lda * 0.5d0
end

subroutine routine_write
 implicit none
 integer :: i,npts
 double precision :: xmax,dx,x
 double precision, allocatable :: dm_hf_array(:), dm_array(:), g0_array(:), on_top_hf(:),on_top_psi(:),on_top_g0(:),mu_hf(:)
 double precision, allocatable :: on_top_exact(:),grad_dm(:)
 double precision :: sqpi
 sqpi = dsqrt(dacos(-1.d0))
 npts = 500
 xmax = 5.d0
 allocate( dm_hf_array(npts), dm_array(npts), g0_array(npts), on_top_psi(npts), on_top_hf(npts), on_top_g0(npts), mu_hf(npts), on_top_exact(npts),grad_dm(npts))
 call routine_psi_hf(dm_hf_array,dm_array,g0_array,on_top_psi,on_top_hf,on_top_g0, mu_hf,grad_dm, npts,xmax)
 call routine_exact(on_top_exact,npts,xmax)
 dx = xmax/dble(npts)
 x = 0.d0
 double precision :: accu1, accu2, accu3, mu_lda, mu_exact
 accu1 = 0.d0
 accu2 = 0.d0
 accu3 = 0.d0
 do i = 1, npts
!  mu_lda   =  - 1.d0 / (dlog(on_top_g0(i)/on_top_hf(i)) * sqpi) 
  mu_lda   =  - 1.d0 / (dlog(2.d0 * g0_array(i)) * sqpi) 
!  mu_lda   =  - 1.d0 / (dlog( g0_array(i) * dm_array(i)**2.d0/(2.d0 * on_top_hf(i))) * sqpi) 
  mu_exact =  - 1.d0 / (dlog(on_top_exact(i)/on_top_hf(i)) * sqpi)
  !                          1     2             3             4                5            
  write(33,'(100(F16.10,X))')x,on_top_psi(i),on_top_hf(i), on_top_g0(i), on_top_exact(i), &  
                             ! 6 
                              on_top_hf(i) + 0.5d0 * dm_array(i)**2.d0 * (g0_array(i) - 0.5d0), &
                             grad_dm(i)/dm_array(i)**(4.d0/3.d0), &
                             !                    6 : n2xc exact 
                             (2.d0 * on_top_exact(i) - dm_array(i)**2.d0)/dm_array(i)**2.d0, & 
                             !    7 : 
                             (2.d0 * on_top_g0(i) - dm_array(i)**2.d0)/dm_array(i)**2.d0
!                             ! 8 !
!                               on_top_exact(i) * (1.d0 + 2.d0/(sqpi * mu_hf(i)))        , & 
!                             !  7              8                                        9              10
!                             g0_array(i), 2.d0 * on_top_exact(i)/ dm_array(i)**2.d0, dm_hf_array(i),dm_array(i), &
!                             !     11                                               12
!                             - 1.d0 / (dlog(on_top_exact(i)/on_top_psi(i))* sqpi) , mu_hf(i), &
!                             !     13   14 
!                             mu_exact, mu_lda
  x += dx
  accu1 += dm_hf_array(i) * dx 
  accu2 += dm_hf_array(i) * dx * mu_lda
  accu3 += dm_hf_array(i) * dx * mu_exact
 enddo
 print*,'mu_exact = ',accu3/accu1
 print*,'mu_lda   = ',accu2/accu1
 print*,'accu1    = ',accu1
end

subroutine routine_psi_hf(dm_hf_array,dm_array,g0_array,on_top_psi,on_top_hf,on_top_g0, mu_hf,grad_dm,  npts,xmax)
 implicit none
 double precision, intent(in) :: xmax
 integer, intent(in) :: npts
 double precision, intent(out) :: dm_hf_array(npts), dm_array(npts), g0_array(npts), on_top_psi(npts), on_top_hf(npts), on_top_g0(npts), mu_hf(npts),grad_dm(npts)
 double precision :: g0_UEG_mu_inf,g0,on_top_in_r
 double precision :: n2_hf,rho_a,rho_b,rho_a_hf,rho_b_hf
 double precision :: mos_array(mo_num)
 double precision :: dx,r(3)
 integer :: i,istate,m
 double precision :: sqpi,f_HF_val_ab,two_bod_dens,w_hf
 double precision ::  dm_a_grad(3), dm_b_grad(3), aos_array(ao_num), grad_aos_array(3, ao_num)
 sqpi = dsqrt(dacos(-1.d0))
 istate = 1
 dx = xmax/dble(npts)
 r = 0.d0
 do i = 1, npts
  call give_all_mos_at_r(r,mos_array)
  call dm_dft_alpha_beta_at_r(r,rho_a,rho_b)
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b,  dm_a_grad, dm_b_grad, aos_array, grad_aos_array)
  call give_on_top_in_r_one_state(r,istate,on_top_in_r)
  on_top_psi(i) = on_top_in_r
  rho_a_hf = mos_array(1)**2.d0
  rho_b_hf = mos_array(1)**2.d0
  dm_hf_array(i) = rho_a_hf + rho_b_hf
  on_top_hf(i) = rho_a_hf * rho_b_hf
  dm_array(i) = rho_a + rho_b
  g0 = g0_UEG_mu_inf(rho_a,rho_b)
  g0_array(i)  = g0 
  grad_dm(i) = 0.d0
  do m = 1, 3
   grad_dm(i) += dm_a_grad(m)**2.d0 + dm_b_grad(m)**2.d0 + dm_a_grad(m) * dm_b_grad(m)
  enddo
!  on_top_g0(i) = g0_array(i) * 0.5d0 * dm_hf_array(i)**2.d0
  call f_HF_valence_ab(r,r,f_HF_val_ab,two_bod_dens)
  w_hf = f_HF_val_ab/two_bod_dens
  mu_hf(i)  = w_hf * sqpi * 0.5d0
  on_top_g0(i) = g0_array(i) * 0.5d0 * dm_array(i)**2.d0
!  on_top_g0(i) = g0_array(i) * 1.0d0 * dm_array(i)**2.d0
  r(1) += dx
 enddo
end


subroutine routine_exact(on_top_exact,npts,xmax)
 implicit none
 double precision, intent(in) :: xmax
 integer, intent(in) :: npts
 double precision, intent(out) :: on_top_exact(npts)
 do i = 1, N_det
  psi_coef(i,1) = reigvec_trans(i,1)/dsqrt(reigvec_trans_norm(1))
 enddo
 touch psi_coef
 double precision :: psi,inv_norm,mu,r12,full_jastrow_mu
 mu = mu_erf 
 print*,'<Phi | e^{2 j(r12,mu)} | Phi> = ',norm_n2_jastrow_cst_mu(1)
 inv_norm = 1.d0/dsqrt(norm_n2_jastrow_cst_mu(1))
 double precision :: dx,r(3)
 integer :: i,istate
 istate = 1
 dx = xmax/dble(npts)
 r = 0.d0
 do i = 1, npts
  call get_two_e_psi_at_r1r2(r,r,psi)
  on_top_exact(i) = (psi * full_jastrow_mu(mu,0.d0) * inv_norm)**2.d0
  r(1) += dx
 enddo
end

 BEGIN_PROVIDER [double precision, average_mu_lda      ]
 implicit none
 integer :: ipoint,i
 double precision :: sqpi
 double precision :: weight, rho_a_hf, rho_b_hf, g0
 double precision :: g0_UEG_mu_inf
 average_mu_lda       = 0.d0
 sqpi = dsqrt(dacos(-1.d0))
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  rho_a_hf = 0.d0
  do i = 1, elec_alpha_num 
   rho_a_hf += mos_in_r_array(i,ipoint)*mos_in_r_array(i,ipoint)
  enddo
  rho_b_hf = 0.d0
  do i = 1, elec_beta_num 
   rho_b_hf += mos_in_r_array(i,ipoint)*mos_in_r_array(i,ipoint)
  enddo
  g0 = g0_UEG_mu_inf(rho_a_hf,rho_b_hf)
  average_mu_lda +=  - 1.d0 / (dlog(2.d0 * g0) * sqpi) * weight * (rho_b_hf+rho_a_hf)
 enddo 
 average_mu_lda = average_mu_lda / dble(elec_alpha_num + elec_beta_num)

END_PROVIDER 
