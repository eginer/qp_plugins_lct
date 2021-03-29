program fcidump_tc_h_3_body_sym
 implicit none
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 call fcidump_3_tc
 call fcidump_2_tc
end

subroutine fcidump_3_tc
 implicit none
 integer :: i,j,k,l,m,n
 integer :: ii,jj,kk,ll,mm,nn
 double precision :: integral 
 character*(128) :: output_physicist
 integer :: i_unit_output_physicist
 integer :: getUnitAndOpen
 output_physicist =trim(ezfio_filename)//'/FCIDUMP_3_body_tc_sym'
 i_unit_output_physicist = getUnitAndOpen(output_physicist,'w')
  do nn = 1, n_act_orb
   n = list_act(nn)
   do ll = 1, n_act_orb
    l = list_act(ll)
    do kk = 1, n_act_orb
     k = list_act(kk)
     do mm = n, n_act_orb
      m = list_act(mm)
      do jj = l, n_act_orb
       j = list_act(jj)
       do ii = k, n_act_orb
        i = list_act(ii)
        !                          1 2 3 1 2 3
        !                         <i j m|k l n>
        !                         (ik|jl|mn)
!        integral = three_body_ints(i,j,m,k,l,n)
        call give_integrals_3_body(i,j,m,k,l,n,integral)
        integral = -integral * 1.d0/3.d0 !!!! For NECI convention 
        if(dabs(integral).lt.1.d-12)cycle
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, i, j, m, k, l, n
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, k,j,m,i,l,n
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, k,l,m,i,j,n
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, k,j,n,i,l,m
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, k,l,n,i,j,m
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, i,l,m,k,j,n
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, k,l,m,i,j,n
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, i,l,n,k,j,m
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, i,j,n,k,l,m
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, k,j,n,i,l,m
        write(i_unit_output_physicist,'(E20.10,6(I3,X))') integral, i,l,n,k,j,m
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
end

subroutine fcidump_2_tc
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
  PROVIDE mo_two_e_integrals_erf_in_map mo_two_e_eff_dr12_pot_array_physicist mo_two_e_integrals_eff_pot_in_map

  double precision :: get_mo_two_e_integral_erf, integral,mo_two_e_integral_eff_pot

  do l=1,n_act_orb
   l1 = list_act(l)
   do k=1,n_act_orb
    k1 = list_act(k)
    do j=1,n_act_orb
     j1 = list_act(j)
     do i=1,n_act_orb
      i1 = list_act(i)
          integral = get_mo_two_e_integral_erf(i1,j1,k1,l1,mo_integrals_erf_map)
          integral += mo_two_e_integral_eff_pot(i1,j1,k1,l1)
          integral += mo_two_e_eff_dr12_pot_array_physicist(i1,j1,k1,l1)
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
