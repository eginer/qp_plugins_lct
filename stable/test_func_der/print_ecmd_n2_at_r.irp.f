program print_ecmd_n2_at_r
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
 call routine_print_intermediaire
end program

subroutine routine_print_intermediaire
 implicit none
 double precision :: r(3),ecmd_pben2_at_r, decdrho_at_r, decdrho2_at_r,decdrho2_at_r_fact, rho2, rho2_extrap,rho,mu,ecPBE

  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.su-ot_print'
  i_unit_output = getUnitAndOpen(output,'w')

 double precision :: xmax, dx
 integer :: nx, ipoint

 r(1) = nucl_coord(1, 1) ! # of the nucleus in the .xyz file
 r(2) = nucl_coord(1, 2)
 r(3) = nucl_coord(1, 3)
 print*,'Nucl coord 1 '
 print*,nucl_coord(1,:)
 print*,'Nucl coord 2 '
 print*,nucl_coord(2,:)
 
 nx = 2000
 xmax = 15.d0
 dx = xmax/dble(nx)
 r(3) += -xmax*0.5d0
 
 write(i_unit_output,*)'#r(3)  ecmd_n2   decdrho  decdrho2  decdrho2*dn2_extrap/dn2   rho2 rho2_extrap  rho  mu'
 do ipoint=0, nx
  call energy_xc_pben2_test_at_r (r,ecmd_pben2_at_r, decdrho_at_r, decdrho2_at_r,decdrho2_at_r_fact, rho2, rho2_extrap,rho,mu,ecPBE)
 write(i_unit_output,'(100(F16.10,X))') r(3),ecmd_pben2_at_r, decdrho_at_r, decdrho2_at_r,decdrho2_at_r_fact,rho2,rho2_extrap,rho,mu,ecPBE
  r(3) += dx
 enddo
end program


subroutine energy_xc_pben2_test_at_r (r,ecmd_pben2_at_r, decdrho_at_r, decdrho2_at_r, decdrho2_at_r_fact, rho2, rho2_extrap,rho,mu,ecPBE)


 implicit none
 BEGIN_DOC
! func_ecmd_utils/rout_pbe_ueg.irp.f
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out):: ecmd_pben2_at_r, decdrho_at_r, decdrho2_at_r, decdrho2_at_r_fact 
 double precision, intent(out):: rho2_extrap,rho2,rho(N_states),mu,ecPBE
 double precision :: rho_a(N_states),rho_b(N_states),grad_rho_a(3, N_states),grad_rho_b(3, N_states)
 double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, aos_array(ao_num), grad_aos_array(3,ao_num)
 double precision :: decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b,ec, decdgrad_rho_2,decdrho2, decdrho
 integer :: i, m, istate
 double precision :: f_psi, dn2_extrap_dn2 
 double precision :: on_top_extrap, mu_correction_of_on_top

   istate = 1
   call give_mu_of_r_cas(r,istate,mu,f_psi,rho2)
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array) 
   rho = rho_a + rho_b 
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m,istate) * grad_rho_a(m,istate)
    grad_rho_b_2 += grad_rho_b(m,istate) * grad_rho_b(m,istate)
    grad_rho_a_b += grad_rho_a(m,istate) * grad_rho_b(m,istate)
   enddo
   call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_2,grad_rho_o_c)
   double precision :: rho_o,rho_c,grad_rho_o_2,grad_rho_2,grad_rho_o_c
   double precision :: decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_o
   call ec_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_2,grad_rho_o_2,grad_rho_o_c,ecPBE,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_o)


!!!!!Dans le main
!  We take the extrapolated on-top pair density (Eq. 29)
!   on_top = total_cas_on_top_density(1,1) !! C'EST PAS LE BON MU ICI
!  Multiplied by 2 because of difference of normalizations between the on_top of QP2 and that of JCP, 150, 084103 1-10 (2019)
   rho2_extrap = 2.d0 * mu_correction_of_on_top(mu,rho2)
   rho2 = 2.d0 * rho2
   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2_extrap,ec,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

   ecmd_pben2_at_r      = ec
   decdrho_at_r         = decdrho
   decdrho2_at_r        = decdrho2
   dn2_extrap_dn2       = 1.d0/(1.d0 + 2.d0/(dacos(-1.d0)**(-0.5d0) * mu)) ! we miss the mu
   decdrho2_at_r_fact   = decdrho2 * dn2_extrap_dn2

end

