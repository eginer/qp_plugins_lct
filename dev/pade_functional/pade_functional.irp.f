program pade_functional
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .true.
  touch read_wf
  io_mo_one_e_integrals = "None"
  touch io_mo_one_e_integrals
  io_mo_two_e_integrals = "None"
  touch io_mo_two_e_integrals
  io_ao_two_e_integrals = "None"
  touch io_ao_two_e_integrals
 
  io_mo_integrals_e_n = "None"
  touch io_mo_integrals_e_n
  io_mo_integrals_kinetic = "None"
  touch io_mo_integrals_kinetic 
  io_ao_integrals_e_n = "None"
  touch io_ao_integrals_e_n 
  io_ao_integrals_kinetic = "None"
  touch io_ao_integrals_kinetic 
! call pote
  call pote_hartree
end

 subroutine pote
 implicit none
 print*,'****************************'
 print*,'psi_energy       =',psi_energy_two_e
 print*,'psi_energy_erf   =',psi_energy_erf
 print*,'****************************'


 print*,'****************************'
 print*,'E_Hx,md^sr_exact =',psi_energy_two_e-psi_energy_erf
!print*,'E_x,md^lr_exact  =',psi_energy_erf-long_range_hartree
 print*,'E_H^sr           =',short_range_Hartree
 print*,'****************************'
 print*,'mu_erf         =',mu_erf
 print*,'e_hx_sr_exp0   =',e_hx_sr_exp0
 print*,'e_hx_sr_exp1   =',e_hx_sr_exp1
 print*,'e_hx_sr_exp2   =',e_hx_sr_exp2
 print*,'e_hx_sr_approx =',e_hx_sr_tot 
 print*,'****************************'
 print*,'Pade 1         =',pade_1
 print*,'Pade 2         =',pade_2
 print*,'Pade 3         =',pade_3
 print*,'****************************'
 print*,' '
 print*,' '
 print*,'****************************'
 print*,'e_hx_local     =',e_hx_local 
 print*,'****************************'
 print*,'Pade 1 local   =',e_hx_pade_local_1
 print*,'Pade 2 local   =',e_hx_pade_local_2
 print*,'Pade 3 local   =',e_hx_pade_local_3
 print*,'****************************'

 end

 subroutine pote_hartree
 implicit none
 print*,'****************************'
 print*,'E_H^sr           =',short_range_Hartree
 print*,'****************************'
 print*,'mu_erf         =',mu_erf
 print*,'e_h_sr_exp0   =',e_h_sr_exp0
 print*,'e_h_sr_exp1   =',e_h_sr_exp1
 print*,'e_h_sr_exp2   =',e_h_sr_exp2
 print*,'e_h_sr_approx =',e_h_sr_tot 
 print*,'****************************'
!print*,'Pade 1         =',pade_1
!print*,'Pade 2         =',pade_2
!print*,'Pade 3         =',pade_3
!print*,'****************************'
!print*,' '
!print*,' '
!print*,'****************************'
!print*,'e_hx_local     =',e_hx_local 
!print*,'****************************'
!print*,'Pade 1 local   =',e_hx_pade_local_1
!print*,'Pade 2 local   =',e_hx_pade_local_2
!print*,'Pade 3 local   =',e_hx_pade_local_3
!print*,'****************************'

 end
