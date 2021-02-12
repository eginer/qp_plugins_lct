program print_ecmd_ueg_at_r
 implicit none
 BEGIN_DOC
 !
 END_DOC

 read_wf = .true.
 touch read_wf
! basis_cor_func = "su_pbe_ot"
! touch basis_cor_func

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

 double precision :: r(3),ecmd_pbeueg_at_r, decdrho_at_r
 double precision :: xmax, dx
 integer :: nx, ipoint
 print *, '#r(3)  ecmd_ueg   decdrho'


 r(1) = nucl_coord(1, 1) ! # of the nucleus in the .xyz file
 r(2) = nucl_coord(1, 2)
 r(3) = nucl_coord(1, 3)
 
 nx = 500
 xmax = 7.d0
 dx = xmax/dble(nx)
 r(3) += -xmax*0.5d0

 do ipoint=0, nx
  call energy_xc_pbeueg_test_at_r (r,ecmd_pbeueg_at_r, decdrho_at_r)
  print*, r(3),ecmd_pbeueg_at_r, decdrho_at_r


  r(3) += dx
 enddo
end program


subroutine energy_xc_pbeueg_test_at_r (r,ecmd_pbeueg_at_r, decdrho_at_r)
 implicit none
 BEGIN_DOC
! func_ecmd_utils/rout_pbe_ueg.irp.f
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out):: ecmd_pbeueg_at_r, decdrho_at_r
 double precision :: rho_a(N_states),rho_b(N_states),grad_rho_a(3, N_states),grad_rho_b(3, N_states)
 double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, aos_array(ao_num), grad_aos_array(3,ao_num)
 double precision :: dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,ex
 double precision :: decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b,ec
 integer :: i, m, istate

   istate = 1
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array) 
   
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m,istate) * grad_rho_a(m,istate)
    grad_rho_b_2 += grad_rho_b(m,istate) * grad_rho_b(m,istate)
    grad_rho_a_b += grad_rho_a(m,istate) * grad_rho_b(m,istate)
   enddo

   call exc_dexc_md_sr_PBE(mu_erf_dft,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

   ecmd_pbeueg_at_r = ec
   decdrho_at_r    = decdrho_a + decdrho_b

end

