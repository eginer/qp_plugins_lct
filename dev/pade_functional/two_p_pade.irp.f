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
           
                                 
 call print_E_H_sr_ex_ground_state

 call tracer_d_E_H
 call tracer_d_E_H_2
 call tracer_two_p_pade_h_4
end


 subroutine tracer_d_E_H
 implicit none
 integer :: n_mu
 double precision :: mu_loc,two_pade_dh_3(n_states),two_pade_dh_4(n_states),two_pade_dh_5(n_states),two_pade_dh_6(n_states),two_pade_dh_7(n_states)

 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 print*,trim(ezfio_filename)

 output=trim(ezfio_filename)//'.pade_der_E_H'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')

 do n_mu = 0,1000
  mu_loc=0.0001d0+dble(n_mu)*0.1d0  

  call give_two_p_pade_dh_3(mu_loc,two_pade_dh_3)
  call give_two_p_pade_dh_4(mu_loc,two_pade_dh_4)
  call give_two_p_pade_dh_5(mu_loc,two_pade_dh_5)
  call give_two_p_pade_dh_6(mu_loc,two_pade_dh_6)
 !call give_two_p_pade_dh_7(mu,two_pade_dh_7)


  write(i_unit_output,'(10(F16.10,X))')mu_loc,two_pade_dh_3(1),two_pade_dh_4(1),two_pade_dh_5(1),two_pade_dh_6(1)
 enddo


end


 subroutine tracer_d_E_H_2
 implicit none
 integer :: istate,i,n_mu
 double precision :: f_mu_1,f_mu_2,d_EH_num,h2,mu_loc
 
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 print*,trim(ezfio_filename)

 output=trim(ezfio_filename)//'.num_der_E_H'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')

 h2 = 0.000000001d0

 do n_mu = 0,1000
  mu_loc= 0.0001d0 + dble(n_mu)*0.1d0  

  call clear_mo_erf_map
  call clear_ao_erf_map
  mu_erf = mu_loc + h2
  touch mu_erf
  mu_erf_dft = mu_erf 
  touch mu_erf_dft
  print*,'mu_erf/mu_erf_dft                 =',mu_erf,'/',mu_erf_dft
  call donner_E_H_sr(f_mu_1)
 
  call clear_mo_erf_map
  call clear_ao_erf_map
  mu_erf = mu_erf - 2.d0 * h2
  touch mu_erf
  mu_erf_dft = mu_erf 
  touch mu_erf_dft
  print*,'mu_erf/mu_erf_dft                 =',mu_erf,'/',mu_erf_dft
  call donner_E_H_sr(f_mu_2)
  print*,'f(x+h),f(x-h)',f_mu_1,f_mu_2
  d_EH_num = (f_mu_1 - f_mu_2)/(2.d0 * h2) 
 
  write(i_unit_output,'(10(F16.10,X))')mu_loc,d_EH_num
 enddo


end



 subroutine tracer_two_p_pade_h_4
 implicit none
 integer :: istate,i,n,n_mu
 double precision :: a,b,h,mu,mu_c,mu_loc
 double precision :: two_pade_dh_4(n_states)
 double precision :: integral, error


 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 print*,trim(ezfio_filename)

 output=trim(ezfio_filename)//'.2p_pade_4_E_H'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')

 do n_mu = 0,1000
  mu_loc=0.0001d0+dble(n_mu)*0.1d0  


  b=0.99999999d0
  a =sqrt(mu_loc**2.d0/(1.d0+mu_loc**2.d0))
  n= 10000000
  double precision, allocatable :: listx(:), listy(:)
  allocate(listx(n),listy(n))
 
  h = (b - a)/dble(n)
 
  do istate = 1, n_states
 
   do i=1,n
    mu = a+dble(i)*h
    mu_c =sqrt(mu**2.d0/(1-mu**2.d0))
    call give_two_p_pade_dh_4(mu_c,two_pade_dh_4)
    listx(i) = mu
    listy(i) = two_pade_dh_4(istate)*(1.d0-mu**2.d0)**(-3.d0/2.d0)
   enddo
 
 
   call DsplIntegr (listx, listy, n , a, b, integral, error)
  
   two_p_pade_h_4(istate) = -integral 
 
  enddo
  write(i_unit_output,'(10(F16.10,X))')mu_loc,two_p_pade_h_4(1)
  deallocate(listx,listy)
 enddo
 end



 subroutine print_E_H_sr_ex_ground_state   
 implicit none

 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 integer :: n_mu,istate,n
 double precision:: mu_loc,f_mu_1
 print*,trim(ezfio_filename)

 output=trim(ezfio_filename)//'.exact_E_H_sr'
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')

 do istate = 1, n_states
  do n_mu = 0,1000
   mu_loc=0.0001+dble(n_mu)*0.1  
   call clear_mo_erf_map
   call clear_ao_erf_map
   mu_erf = mu_loc 
   touch mu_erf
   mu_erf_dft = mu_erf 
   touch mu_erf_dft
   call donner_E_H_sr(f_mu_1)
   write(i_unit_output,'(10(F16.10,X))')mu_loc,f_mu_1
  enddo
 enddo  

 end


 subroutine donner_E_H_sr(tmp)
 implicit none
 double precision, intent(out) :: tmp 
 tmp=short_range_Hartree(1)
!print*,''
!print*,'sr Hratree dans la subrouteen',tmp
!print*,'mu_erf  = ',mu_erf
!print*,''
 end
