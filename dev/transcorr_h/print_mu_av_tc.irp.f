program print_mu_av_tc
 implicit none
 read_wf = .True.
 touch read_wf
   print*,'average_mu_lda       = ',average_mu_lda
   print*,'average_mu_rs        = ',average_mu_rs 
   print*,'average_mu_rs_c      = ',average_mu_rs_c
   print*,'average_mu_rs_c_lda  = ',average_mu_rs_c_lda
   print*,'average_mu_grad_n    = ',average_mu_grad_n
   call plot_mu_tc
end


 BEGIN_PROVIDER [double precision, average_mu_lda      ]
&BEGIN_PROVIDER [double precision, average_mu_rs       ]
&BEGIN_PROVIDER [double precision, average_mu_rs_c     ]
&BEGIN_PROVIDER [double precision, average_mu_rs_c_lda ]
&BEGIN_PROVIDER [double precision, average_mu_grad_n   ]
 implicit none
 integer :: ipoint,i,m
 double precision :: sqpi
 double precision :: weight, rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs
 sqpi = dsqrt(dacos(-1.d0))
 cst_rs   = (4.d0 * dacos(-1.d0)/3.d0)**(-1.d0/3.d0)
 alpha_rs = 2.d0 * dsqrt((9.d0 * dacos(-1.d0)/4.d0)**(-1.d0/3.d0)) / sqpi
 average_mu_lda     = 0.d0
 average_mu_rs      = 0.d0
 average_mu_rs_c    = 0.d0
 average_mu_grad_n  = 0.d0 
 double precision :: elec_a,elec_b
 elec_a = 0.d0
 elec_b = 0.d0
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  rho_a_hf = 0.d0
  grad_n   = 0.d0
  rho_a_hf = one_e_dm_and_grad_alpha_in_r(4,ipoint,1)
  rho_b_hf = one_e_dm_and_grad_beta_in_r(4,ipoint,1)
  grad_n = one_e_grad_2_dm_alpha_at_r(ipoint,1) + one_e_grad_2_dm_beta_at_r(ipoint,1)
  grad_n += 2.d0 * scal_prod_grad_one_e_dm_ab(ipoint,1)
  rho_hf = rho_a_hf + rho_b_hf
  grad_n = dsqrt(grad_n)
  if(dabs(rho_hf).gt.1.d-20)then
   grad_n = grad_n/(4.d0 * rho_hf)
  else
   grad_n = 0.d0
  endif
  rs = cst_rs * rho_hf**(-1.d0/3.d0)
  g0 = g0_UEG_mu_inf(rho_a_hf,rho_b_hf)

  average_mu_rs     += 1.d0/rs * rho_hf * weight
  average_mu_rs_c   += alpha_rs/dsqrt(rs) * rho_hf * weight
  average_mu_lda    +=  - 1.d0 / (dlog(2.d0 * g0) * sqpi) * weight * rho_hf
  average_mu_grad_n += grad_n * rho_hf * weight
  elec_a += rho_a_hf * weight
  elec_b += rho_b_hf * weight
 enddo 
 print*,'elec_a,elec_b',elec_a,elec_b
 average_mu_lda    = average_mu_lda   / dble(elec_a+ elec_b)
 average_mu_rs     = average_mu_rs    / dble(elec_a+ elec_b)
 average_mu_rs_c   = average_mu_rs_c  / dble(elec_a+ elec_b)
 average_mu_grad_n = average_mu_grad_n/ dble(elec_a+ elec_b)

 average_mu_rs_c_lda = 0.5d0 * (average_mu_lda + average_mu_rs_c)

END_PROVIDER 

subroutine plot_mu_tc
 implicit none
 integer :: i,nx,m
 double precision :: xmin,xmax,dx
 double precision :: r(3)
 double precision :: sqpi
 double precision :: rho_a_hf, rho_b_hf, g0,rho_hf
 double precision :: rs,grad_n,grad_n_a(3), grad_n_b(3)
 double precision :: g0_UEG_mu_inf
 double precision :: cst_rs,alpha_rs,mu_rs, mu_rs_c, mu_lda, mu_grad_n 
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 if (no_core_density)then
  output=trim(ezfio_filename)//'.tc_mu_fc'
 else
  output=trim(ezfio_filename)//'.tc_mu'
 endif
 i_unit_output = getUnitAndOpen(output,'w')

 nx = 5000
 xmin = -10.d0
 xmax = 10.d0
 dx   = (xmax - xmin)/dble(nx)
 r = 0.d0
 r(3) = xmin

 sqpi = dsqrt(dacos(-1.d0))
 cst_rs   = (4.d0 * dacos(-1.d0)/3.d0)**(-1.d0/3.d0)
 alpha_rs = 2.d0 * dsqrt((9.d0 * dacos(-1.d0)/4.d0)**(-1.d0/3.d0)) / sqpi

  write(i_unit_output,*)'#       r(3)              rho_hf          mu_rs           mu_rs_c           mu_lda         mu_grad_n'
 do i = 1, nx
  rho_a_hf = 0.d0
  grad_n   = 0.d0
  call density_and_grad_alpha_beta(r,rho_a_hf,rho_b_hf, grad_n_a, grad_n_b)
  do m = 1, 3
   grad_n += grad_n_a(m)*grad_n_a(m) + grad_n_b(m)*grad_n_b(m) + 2.d0 * grad_n_a(m) * grad_n_b(m)
  enddo
  rho_hf = rho_a_hf + rho_b_hf
  grad_n = dsqrt(grad_n)
  grad_n = grad_n/(4.d0 * rho_hf)
  rs = cst_rs * rho_hf**(-1.d0/3.d0)
  g0 = g0_UEG_mu_inf(rho_a_hf,rho_b_hf)

  mu_rs     = 1.d0/rs 
  mu_rs_c   = alpha_rs/dsqrt(rs) 
  mu_lda    =  - 1.d0 / (dlog(2.d0 * g0) * sqpi)  
  mu_grad_n = grad_n 
  write(i_unit_output,'(100(F16.10,X))')r(3),rho_hf,mu_rs,mu_rs_c,mu_lda,mu_grad_n
  r(3) += dx
 enddo


end
