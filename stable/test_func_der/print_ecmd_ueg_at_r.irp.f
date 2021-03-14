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

 call routine_print_intermediaire

end program

subroutine routine_print_intermediaire
 implicit none
 double precision :: r(3),ecmd_pbeueg_at_r, decdrho_at_r,rho
 double precision :: xmax, dx,xmin
 integer :: nx, ipoint
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.pbe-ueg_print'
  i_unit_output = getUnitAndOpen(output,'w')


 r(1) = nucl_coord(1, 1) ! # of the nucleus in the .xyz file
 r(2) = nucl_coord(1, 2)
 r(3) = nucl_coord(1, 3)
 print*,'Nucl coord 1 '
 print*,nucl_coord(1,:)
 print*,'Nucl coord 2 '
 print*,nucl_coord(2,:)
 
 xmin = -4.d0
 xmax =  5.d0
 nx = 1000
 dx = (xmax - xmin)/dble(nx)
 r(3) = xmin 

   write(i_unit_output,*)'#r(3)  ecmd_ueg   decdrho  rho'
 do ipoint=0, nx
  call energy_xc_pbeueg_test_at_r (r,ecmd_pbeueg_at_r, decdrho_at_r,rho)
  write(i_unit_output,'(100(F16.10,X))') r(3),ecmd_pbeueg_at_r, decdrho_at_r,rho
  r(3) += dx
 enddo
end


subroutine energy_xc_pbeueg_test_at_r (r,ecmd_pbeueg_at_r, decdrho_at_r,rho)
 implicit none
 BEGIN_DOC
! func_ecmd_utils/rout_pbe_ueg.irp.f
 END_DOC
 double precision, intent(in) :: r(3)
 double precision, intent(out):: ecmd_pbeueg_at_r, decdrho_at_r,rho(N_states)
 double precision :: rho_a(N_states),rho_b(N_states),grad_rho_a(3, N_states),grad_rho_b(3, N_states)
 double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, aos_array(ao_num), grad_aos_array(3,ao_num)
 double precision :: dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,ex
 double precision :: decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b,ec
 integer :: i, m, istate

   istate = 1
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
   call exc_dexc_md_sr_PBE(mu_erf_dft,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

   ecmd_pbeueg_at_r = ec
   decdrho_at_r    = decdrho_a + decdrho_b

end

