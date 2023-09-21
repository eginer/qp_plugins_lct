
subroutine fcidump_2_tc_cst_mu
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
  output=trim(ezfio_filename)//'/FCIDUMP_2_body_tc'
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
  PROVIDE mo_non_hermit_term mo_tc_sym_two_e_pot_in_map

  double precision :: integral,get_mo_tc_sym_two_e_pot

  do l=1,n_act_orb
   l1 = list_act(l)
   do k=1,n_act_orb
    k1 = list_act(k)
    do j=1,n_act_orb
     j1 = list_act(j)
     do i=1,n_act_orb
      i1 = list_act(i)
          integral = get_mo_tc_sym_two_e_pot(i1,j1,k1,l1,mo_tc_sym_two_e_pot_map)
          integral += mo_non_hermit_term(i1,j1,k1,l1)
          if (dabs(integral) > mo_integrals_threshold) then
            write(i_unit_output,*) integral, i,k,j,l
          endif
     enddo
    enddo
   enddo
  enddo

  do j=1,n_act_orb
   j1 = list_act(j)
   do i=1,n_act_orb
    i1 = list_act(i)
      integral = mo_one_e_integrals(i1,j1) + core_fock_operator(i1,j1)
      if (dabs(integral) > mo_integrals_threshold) then
        write(i_unit_output,*) integral, i,j,0,0
      endif
   enddo
  enddo
  write(i_unit_output,*) core_energy, 0, 0, 0, 0
end
