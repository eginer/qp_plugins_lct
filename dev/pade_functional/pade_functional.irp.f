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
  io_mo_two_e_integrals_erf = "None"
  touch io_mo_two_e_integrals_erf
  io_ao_two_e_integrals_erf = "None"
  touch io_ao_two_e_integrals_erf
 
  io_mo_integrals_e_n = "None"
  touch io_mo_integrals_e_n
  io_mo_integrals_kinetic = "None"
  touch io_mo_integrals_kinetic 
  io_ao_integrals_e_n = "None"
  touch io_ao_integrals_e_n 
  io_ao_integrals_kinetic = "None"
  touch io_ao_integrals_kinetic 
! call pote
!  call pote_hartree
! call pote_test
! call pote_oredre_infinityyyyyy
! call pote_oredre_infinityyyyyy_pade
! call pote_hartree_infiny
! call pote_hartree_test_exp_mu_0
  call two_points_pade_hartree_test
! call tracer_d_E_H
!call pote_test
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
 print*,'Pade 1         =',pade_1_h
 print*,'Pade 2         =',pade_2_h
 print*,'Pade 3         =',pade_3_h
 print*,'****************************'
 print*,' '
 print*,' '
 print*,'****************************'
 print*,'e_h_local     =',e_h_local 
 print*,'****************************'
 print*,'Pade 1 local   =',e_h_pade_local_1
 print*,'Pade 2 local   =',e_h_pade_local_2
 print*,'Pade 3 local   =',e_h_pade_local_3
 print*,'****************************'

 end

 subroutine pote_hartree_infiny
 implicit none
 print*,'****************************'
 print*,'E_H^sr           =',short_range_Hartree
 print*,'****************************'
 print*,'mu_erf         =',mu_erf
 print*,'e_h_pade_diag_coef   =',e_h_pade_diag_coef(1)
 print*,'****************************'
 end


 subroutine pote_hartree_test_exp_mu_0
 implicit none
 print*,'****************************'
 print*,'E_H^sr                        =',short_range_Hartree
 print*,'****************************'
 print*,'mu_erf                        =',mu_erf
 print*,'regular_range_Hartree         =',regular_range_Hartree(1)
 print*,'sr_hatree_energy_small_mu_tot =',sr_hatree_energy_small_mu_tot(1)
 print*,'Contriub expension            =',sr_hatree_energy_small_mu_tot(1)-regular_range_Hartree(1)
 print*,'****************************'
 print*,'Pade test                     =',regular_range_Hartree(1) + e_h_pade_mu_0_diag_coef(1)
 print*,'Pade contrib                  =',e_h_pade_mu_0_diag_coef(1)
 end


 subroutine pote_oredre_infinityyyyyy
 implicit none

 print*,'****************************'
 print*,'E_Hx,md^sr_exact =',psi_energy_two_e-psi_energy_erf
 print*,'****************************'
 print*,'mu_erf         =',mu_erf
 print*,'e_hx_sr_exp0   =',e_hx_sr_by_order_exp(1,0)
 print*,'e_hx_sr_exp1   =',e_hx_sr_by_order_exp(1,1)
 print*,'e_hx_sr_exp2   =',e_hx_sr_by_order_exp(1,2)
 print*,'e_hx_sr_exp3   =',e_hx_sr_by_order_exp(1,3)
 print*,'e_hx_sr_exp4   =',e_hx_sr_by_order_exp(1,4)
 print*,'e_hx_sr_exp5   =',e_hx_sr_by_order_exp(1,5)
 print*,'e_hx_sr_exp6   =',e_hx_sr_by_order_exp(1,6)
 print*,'****************************'

 end


 subroutine pote_oredre_infinityyyyyy_pade
 implicit none

 print*,'****************************'
 print*,'E_Hx,md^sr_exact =',psi_energy_two_e-psi_energy_erf
 print*,'****************************'
 print*,'mu_erf         =',mu_erf
 print*,'e_hx_pade_diag_coef   =',e_hx_pade_diag_coef(1)
! veeeerrriiiifffffffffffff
!print*,'e_hx_sr_exp0   =',e_hx_sr_exp0
!print*,'poteteeete',e_hx_sr_by_order_exp(1,1) 
!print*,'poteteeete',e_hx_sr_by_order_exp(1,2) 
!print*,'poteteeete',e_hx_sr_by_order_exp(1,3) 
!print*,'poteteeete',e_hx_sr_by_order_exp(1,4) 
 print*,'pade Old shcool 1 =',pade_1
 print*,'pade Old shcool 2 =',pade_2
 print*,'pade Old shcool 2 =',pade_3
 print*,'****************************'
 end

 subroutine pote_test   
 implicit none
 provide simpson_int_test
 end



 subroutine two_points_pade_hartree_test
 implicit none
 double precision :: two_pade_dh_4(n_states),d_EH_num,h,f_mu_1,f_mu_2 
!provide simpson_int_test
 !rovide two_p_pade_h_3
!provide two_point_pade_h_a_4
 print*,'****************************'
 print*,'Regular_range_Hartree         =',regular_range_Hartree(1)
 print*,'E_H^sr                        =',short_range_Hartree
 print*,'****************************'
 print*,'mu_erf                        =',mu_erf
!print*,"two_p_pade_h_3                =",two_p_pade_h_3(1)
 print*,"two_p_pade_h_4                =",two_p_pade_h_4(1)
!print*,"two_p_pade_h_5                =",two_p_pade_h_5(1)
!print*,"two_p_pade_h_6                =",two_p_pade_h_6(1)
!print*,"two_p_pade_h_7                =",two_p_pade_h_7(1)
 print*,'****************************'
!print*,'Abs error pade 3              =',short_range_Hartree-two_p_pade_h_3(1)
 print*,'Abs error pade 4              =',short_range_Hartree-two_p_pade_h_4(1)
!print*,'Abs error pade 5              =',short_range_Hartree-two_p_pade_h_5(1)
!print*,'Abs error pade 6              =',short_range_Hartree-two_p_pade_h_6(1)
!print*,'Abs error pade 7              =',short_range_Hartree-two_p_pade_h_7(1)
 print*,'****************************'
!print*,'Relat error pade 3            =',(short_range_Hartree-two_p_pade_h_3(1))/short_range_Hartree
 print*,'Relat error pade 4            =',(short_range_Hartree-two_p_pade_h_4(1))/short_range_Hartree
!print*,'Relat error pade 5            =',(short_range_Hartree-two_p_pade_h_5(1))/short_range_Hartree
!print*,'Relat error pade 6            =',(short_range_Hartree-two_p_pade_h_6(1))/short_range_Hartree
!print*,'Relat error pade 7            =',(short_range_Hartree-two_p_pade_h_7(1))/short_range_Hartree
 print*,'****************************'
!call give_two_p_pade_dh_4(mu_erf,two_pade_dh_4)
!print*,'dE_H_2p_pade_4            =',two_pade_dh_4(1)

!!!!!!!!!!!!!derive numeeeeriiiiqqquuuuuueeeeeee!!!
 print*,'****************************'
 call give_two_p_pade_dh_4(mu_erf,two_pade_dh_4)


 h = 0.000000001
 call clear_mo_erf_map
 call clear_ao_erf_map
 mu_erf = mu_erf + h
 touch mu_erf
 mu_erf_dft = mu_erf 
 touch mu_erf_dft
!print*,'mu_erf/mu_erf_dft                 =',mu_erf,'/',mu_erf_dft
!integer :: get_ao_erf_map_size
!print*,'size AO map = ',get_ao_erf_map_size(ao_integrals_erf_map)
 call donner_E_H_sr(f_mu_1)

 call clear_mo_erf_map
 call clear_ao_erf_map
 mu_erf = mu_erf - 2.d0 * h
 touch mu_erf
 mu_erf_dft = mu_erf 
 touch mu_erf_dft
!print*,'mu_erf/mu_erf_dft                 =',mu_erf,'/',mu_erf_dft
!integer :: get_ao_erf_map_size
!print*,'size AO map = ',get_ao_erf_map_size(ao_integrals_erf_map)
 call donner_E_H_sr(f_mu_2)
 print*,'f(x+h),f(x-h)',f_mu_1,f_mu_2
 d_EH_num = (f_mu_1 - f_mu_2)/(2.d0 * h) 

 print*,'d_EH_num     =',two_pade_dh_4(1) 
 print*,'d_EH pade 4  =',d_EH_num
 print*,'Abs error d_Eh           =',d_EH_num-two_pade_dh_4(1)
 print*,'Relat error d_Eh         =',(d_EH_num-two_pade_dh_4(1))/d_EH_num
 print*,'****************************'

 end


 subroutine donner_E_H_sr(tmp)
 implicit none
 double precision, intent(out) :: tmp 
 tmp=short_range_Hartree(1)
 print*,''
 print*,'sr Hratree dans la subrouteen',tmp
 print*,'mu_erf  = ',mu_erf
 print*,''
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
