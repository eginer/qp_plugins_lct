program fcidump_erf
  implicit none
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

 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 
 call routine

end
subroutine routine
  implicit none
  BEGIN_DOC
! Produce a regular `FCIDUMP` file from the |MOs| stored in the |EZFIO|
! directory.
!
! To specify an active space, the class of the |MOs| have to set in the
! |EZFIO| directory (see :ref:`qp_set_mo_class`).
!
! The :ref:`fcidump` program supports 3 types of |MO| classes :
!
! * the *core* orbitals which are always doubly occupied in the
!   calculation
!
! * the *deleted* orbitals that are never occupied in the calculation
!
! * the *active* orbitals that are occupied with a varying number of
!   electrons
!
  END_DOC
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.FCIDUMP_erf'
  i_unit_output = getUnitAndOpen(output,'w')

  integer :: i,j,k,l
  integer :: i1,j1,k1,l1
  integer :: i2,j2,k2,l2
  integer*8 :: m
  character*(2), allocatable :: A(:)

  write(i_unit_output,*) '&FCI NORB=', n_act_orb, ', NELEC=', elec_num-n_core_orb*2, &
   ', MS2=', (elec_alpha_num-elec_beta_num), ','
  allocate (A(n_act_orb))
  A = '1,'
  write(i_unit_output,*) 'ORBSYM=', (A(i), i=1,n_act_orb)
  write(i_unit_output,*) 'ISYM=0,'
  write(i_unit_output,*) '/'
  deallocate(A)

  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  integer(cache_map_size_kind)   :: n_elements, n_elements_max
  PROVIDE mo_two_e_integrals_in_map

  double precision :: get_mo_two_e_integral_erf, integral

  do l=1,n_act_orb
   l1 = list_act(l)
   do k=1,n_act_orb
    k1 = list_act(k)
    do j=l,n_act_orb
     j1 = list_act(j)
     do i=k,n_act_orb
      i1 = list_act(i)
       if (i1>=j1) then
          integral = get_mo_two_e_integral_erf(i1,j1,k1,l1,mo_integrals_erf_map)
          if (dabs(integral) > mo_integrals_threshold) then
            write(i_unit_output,*) integral, i,k,j,l
          endif
       end if
     enddo
    enddo
   enddo
  enddo

  do j=1,n_act_orb
   j1 = list_act(j)
   do i=j,n_act_orb
    i1 = list_act(i)
      integral = effective_one_e_potential(i1,j1,1) + core_fock_operator_erf(i1,j1) 
      if (dabs(integral) > mo_integrals_threshold) then
        write(i_unit_output,*) integral, i,j,0,0
      endif
   enddo
  enddo
 double precision :: core_e
 core_e = core_energy_erf
 do i = 1, n_core_orb
   j = list_core(i)
   core_e += 2.d0 * effective_one_e_potential(j,j,1)
 enddo
  write(i_unit_output,*) core_e, 0, 0, 0, 0
end
