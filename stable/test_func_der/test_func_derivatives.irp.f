program test_energy_integral_lda
 implicit none
 BEGIN_DOC
 !
 END_DOC

 read_wf = .true.
 touch read_wf
 basis_cor_func = "su_pbe_ot"
 touch basis_cor_func

 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 
 call routine 

   print*, '___________LDA___________'
   print*, "ex_lda [n]                        =", ex_lda_at_n
   print*, "ex_lda [n + delta_n]              =", ex_lda_at_n_plus_delta_n
   print*, "ex_lda [n + delta_n] - ex_lda [n] =", ex_lda_at_n_plus_delta_n - ex_lda_at_n
   print*, "n variation integral              :", int_vx_lda_at_n 
   print*, "------> Relative error            :", (int_vx_lda_at_n - ex_lda_at_n_plus_delta_n + ex_lda_at_n) / (ex_lda_at_n_plus_delta_n - ex_lda_at_n)

   print*, '___________PBE___________'
   print*, "ex_pbe [n] =", ex_pbe_at_n
   print*, "ex_pbe [n + delta_n]              =", ex_pbe_at_n_plus_delta_n
   print*, "ex_pbe [n + delta_n] - ex_pbe [n] =", ex_pbe_at_n_plus_delta_n - ex_pbe_at_n
   print*, "n variation integral              :", int_vx_pbe_at_n 
   print*, "------> Relative error            :", (int_vx_pbe_at_n - ex_pbe_at_n_plus_delta_n + ex_pbe_at_n) / (ex_pbe_at_n_plus_delta_n - ex_pbe_at_n)
   print*, "ec_pbe [n]                        =", ec_pbe_at_n
   print*, "ec_pbe [n + delta_n]              =", ec_pbe_at_n_plus_delta_n
   print*, "ec_pbe [n + delta_n] - ex_pbe [n] =", ec_pbe_at_n_plus_delta_n - ec_pbe_at_n
   print*, "n variation integral              :", int_vc_pbe_at_n 
   print*, "------> Relative error            :", (int_vc_pbe_at_n - ec_pbe_at_n_plus_delta_n + ec_pbe_at_n) / (ec_pbe_at_n_plus_delta_n - ec_pbe_at_n)
! print*, '___________PBEUEG___________'
! print*, "ex_pbeUEG [n, grad_n]                                                  =", ex_pbeUEG_at_n
! print*, "ex_pbeUEG [n + delta_n, grad_n + delta_grad_n]                         =", ex_pbeUEG_at_n_plus_delta_n
! print*, "ex_pbeUEG [n + delta_n, grad_n + delta_grad_n] - ex_pbeUEG [n, grad_n] =", ex_pbeUEG_at_n_plus_delta_n - ex_pbeUEG_at_n
! print*, "n and grad_n variations integral                                       :", int_vx_pbeUEG_at_n 
! print*, "ec_pbeUEG [n, grad_n]                                                  =", ec_pbeUEG_at_n
! print*, "ec_pbeUEG [n + delta_n, grad_n + delta_grad_n]                         =", ec_pbeUEG_at_n_plus_delta_n
! print*, "ec_pbeUEG [n + delta_n, grad_n + delta_grad_n] - ex_pbeUEG [n, grad_n] =", ec_pbeUEG_at_n_plus_delta_n - ec_pbeUEG_at_n
! print*, "n and grad_n variations integral                                       :", int_vc_pbeUEG_at_n 

   print*, '___________PBEn2___________'
   print*, "ec_pben2 [n,n2]                                         =", ec_pben2_at_n_n2
   print*, "ec_pben2 [n + delta_n, n2 + delta_n2]                   =", ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2
   print*, "ec_pben2 [n + delta_n, n2]                              =", ec_pben2_at_n_plus_delta_n
   print*, "ec_pben2 [n, n2 + delta_n2]                             =", ec_pben2_at_n2_plus_delta_n2
   print*, "ec_pben2 [n + delta_n, n2 + delta_n2] - ec_pben2 [n,n2] =", ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2 - ec_pben2_at_n_n2
   print*, "n and n2 variation integral                             :", int_vc_pben2_one_e_at_n_n2 + int_vc_pben2_two_e_at_n_n2
   print*, "--------> Relative error                                :", (int_vc_pben2_one_e_at_n_n2 + int_vc_pben2_two_e_at_n_n2 - ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2 + ec_pben2_at_n_n2) / (ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2 - ec_pben2_at_n_n2)
   print*, "ec_pben2 [n + delta_n, n2] - ec_pben2 [n,n2]            =", ec_pben2_at_n_plus_delta_n - ec_pben2_at_n_n2
   print*, " only n variation integral                              :", int_vc_pben2_one_e_at_n_n2
   print*, "--------> Relative error                                :", (int_vc_pben2_one_e_at_n_n2 - ec_pben2_at_n_plus_delta_n + ec_pben2_at_n_n2)/(ec_pben2_at_n_plus_delta_n - ec_pben2_at_n_n2)
   print*, "ec_pben2 [n, n2 + delta_n2] - ec_pben2 [n,n2]           =", ec_pben2_at_n2_plus_delta_n2 - ec_pben2_at_n_n2
   print*, " only n2 variation integral                             :", int_vc_pben2_two_e_at_n_n2 
   print*, "--------> Relative error                                :", (int_vc_pben2_two_e_at_n_n2 - ec_pben2_at_n2_plus_delta_n2 + ec_pben2_at_n_n2)/(ec_pben2_at_n2_plus_delta_n2 - ec_pben2_at_n_n2)

end

subroutine routine
 implicit none
 include 'constants.include.F'
  integer :: nx, p, m
  double precision :: xmax, dx, r(3), mu, rho_2, f_psi, on_top_extrap, factor, mu_correction_of_on_top
  double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: aos_array(ao_num), grad_aos_array(ao_num,3)
  double precision :: ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2
  double precision :: mos_array(mo_num)

  nx = 500
  xmax = 2.d0
  dx = xmax/dble(nx)
  r(:) = nucl_coord_transp(:,1)
  r(3) += - xmax * 0.5d0
  write(44,*) '#1=r(3) 2=decdrho2 3=factor 4=product'
  do p=1, nx
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo
   call give_all_mos_at_r(r,mos_array)
   
   call give_mu_of_r_cas(r,1,mu,f_psi,rho_2)
   rho_2 = rho_2*2.d0
   on_top_extrap =  mu_correction_of_on_top(mu,rho_2) 
   
   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho_2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

   factor = 1.d0/(1.d0 + 2.d0/(sqpi*mu))
   !                            1    2            3           4         5           6      7     8            9                         10
   write(44,'(100(F16.10,X))') r(3), decdrho2, factor, decdrho2*factor,rho_2,on_top_extrap,mu,rho_a+rho_b, 2.d0 * mos_array(1)**2.d0,mos_array(1)**4.d0
   r(3) += dx
  enddo
end program

