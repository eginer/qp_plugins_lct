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

  call print_exp_big_mu_ground_state
  call print_exp_pade_mu_ground_state
! call print_E_H_sr_ex_ground_state
end

 subroutine print_exp_big_mu_ground_state   
 implicit none

 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 integer :: n_mu,istate,n
 double precision:: mu_loc,e_h_sr_tot_exp_loc 
 print*,trim(ezfio_filename)
 character*(128) :: filename
 write (filename, "(I1)")order_derivative_ontop

 output=trim(ezfio_filename)//'.big_mu_exp_order_'//trim(filename)
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')

 do istate = 1, n_states
  do n_mu = 0,1000
   mu_loc=0.0001+dble(n_mu)*0.1  
   e_h_sr_tot_exp_loc = 0.d0
   do n = 0, order_derivative_ontop      
    e_h_sr_tot_exp_loc += coeff_H_exp(n,istate)/(mu_loc**(2.d0*dble(n)+2.d0))
   enddo
   write(i_unit_output,'(10(F16.10,X))')mu_loc,e_h_sr_tot_exp_loc
  enddo
 enddo  


 end


 subroutine print_exp_pade_mu_ground_state   
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: n,m,istate,INFO,i,j
 n=(order_derivative_ontop+1)/2  
 m=order_derivative_ontop+1-n
 print*,"n =",n
 print*,"m =",m

 double precision :: denom,nume,x_loc,e_h_pade_diag_coef_loc
 double precision, allocatable :: matrix_alpha(:,:)
 double precision, allocatable :: vect_b_temp(:),vect_a(:),vect_b(:),c_pade(:)
 integer, allocatable :: IPIV(:)
 allocate(matrix_alpha(m,m),vect_b_temp(m),IPIV(m),vect_a(0:n),vect_b(0:m+1),c_pade(0:n+m+1))

 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 provide ezfio_filename 
 integer :: n_mu
 double precision:: mu_loc
 print*,trim(ezfio_filename)
 character*(128) :: filename
 write (filename, "(I1)")order_derivative_ontop

 output=trim(ezfio_filename)//'.big_mu_pade_order_'//trim(filename)
 output=trim(output)
 print*,'output = ',trim(output)
 i_unit_output = getUnitAndOpen(output,'w')

 do istate = 1, n_states
 matrix_alpha= 0.d0
 vect_b_temp=0.d0

  c_pade = 0.d0
  do i = 1,order_derivative_ontop+1
   c_pade(i)=coeff_H_exp(i-1,istate)
  enddo 


  do j =1,m
   vect_b_temp(j) = - c_pade(n+j)
   do i = 1,m
   !matrix_alpha(i,j)= 0.d0
    matrix_alpha(i,j)=c_pade(n+i-j)
   !print*,i,j,matrix_alpha(i,j)
   !print*,j,vect_b_temp(j)
   enddo
  enddo
 
  call dgesv(m,1,matrix_alpha,m,IPIV,vect_b_temp,m,INFO) 
 
  print*,'***********************'  
  print*,'INFO=',INFO  
  print*,'***********************'  

  vect_b = 0.d0
  vect_b(0) = 1.d0
  do i=1,m
   vect_b(i)=vect_b_temp(i)
  enddo

  vect_a(0) = 0.d0
  do i = 1,n
   do j = 0,i
    vect_a(i) += c_pade(i-j)*vect_b(j)
   enddo
  enddo

  do n_mu = 0,1000
   mu_loc=0.0001+dble(n_mu)*0.1  

   x_loc=1/(mu_loc**2.d0) 
 
   nume = 0.d0
   do i= 0,n
    nume += vect_a(i)*x_loc**(dble(i))
   enddo
   
   denom = 0.d0
   do i= 0,m
    denom += vect_b(i)*x_loc**(dble(i))
   enddo
  
   e_h_pade_diag_coef_loc = nume/denom

   write(i_unit_output,'(10(F16.10,X))')mu_loc,e_h_pade_diag_coef_loc
  enddo
 
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
