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

  call pote_oredre_Hartree
end


 subroutine pote_oredre_Hartree
 implicit none

 print*,'****************************'
 print*,'E_H^sr           =',short_range_Hartree
 print*,'****************************'
 print*,'mu_erf         =',mu_erf
 print*,'e_h_sr_exp0   =',e_h_sr_by_order_exp(1,0)
 print*,'e_h_sr_exp1   =',e_h_sr_by_order_exp(1,1)
 print*,'e_h_sr_exp2   =',e_h_sr_by_order_exp(1,2)
 print*,'e_h_sr_exp3   =',e_h_sr_by_order_exp(1,3)
 print*,'e_h_sr_exp4   =',e_h_sr_by_order_exp(1,4)
 print*,'e_h_sr_exp5   =',e_h_sr_by_order_exp(1,5)
 print*,'e_h_sr_exp6   =',e_h_sr_by_order_exp(1,6)
 print*,'****************************'

 end

