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

 !call pote_oredre_Hartree
  call tracer_d_E_H 
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

 subroutine tracer_d_E_H
 implicit none
 integer :: istate,i,n
 double precision :: mu,a,b,h,two_pade_dh_3(n_states),two_pade_dh_4(n_states),two_pade_dh_5(n_states),two_pade_dh_6(n_states),two_pade_dh_7(n_states)
 b=40.d0
 a = 0.001d0
 n= 100000
 h = (b - a)/dble(n)
 
 mu = a
 do i=1,n
  mu += h 
  call give_two_p_pade_dh_3(mu,two_pade_dh_3)
  call give_two_p_pade_dh_4(mu,two_pade_dh_4)
  call give_two_p_pade_dh_5(mu,two_pade_dh_5)
  call give_two_p_pade_dh_6(mu,two_pade_dh_6)
  call give_two_p_pade_dh_7(mu,two_pade_dh_7)
  write(33,'(10(F16.10,X))')mu,two_pade_dh_3(1),two_pade_dh_4(1),two_pade_dh_5(1),two_pade_dh_6(1),two_pade_dh_7(1)
 enddo


end
