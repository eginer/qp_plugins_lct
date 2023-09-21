
subroutine FOBOCI_lmct_mlct_old_thr(iter)
  use bitmasks
 implicit none
 integer, intent(in) :: iter
 integer :: i,j,k,l
 integer(bit_kind),allocatable :: unpaired_bitmask(:,:)
 integer, allocatable :: occ(:,:)
 integer       :: n_occ_alpha, n_occ_beta
 double precision :: norm_tmp(N_states),norm_total(N_states)
 logical :: test_sym
 double precision :: thr,hij
 double precision, allocatable :: dressing_matrix(:,:)
 logical :: verbose,is_ok,is_ok_perturbative
 integer(bit_kind) , allocatable :: zero_bitmask(:,:)
 integer(bit_kind) , allocatable :: psi_singles(:,:,:)
 logical :: lmct,exit_loop
 double precision, allocatable :: psi_singles_coef(:,:)
 integer :: i_hole_osoci, i_particl_osoci
 double precision :: e_pt2,  hkl, accu
 verbose = .True.
 thr = 1.d-12
 allocate(unpaired_bitmask(N_int,2))
 allocate (occ(N_int*bit_kind_size,2))
 do i = 1, N_int
  unpaired_bitmask(i,1) = unpaired_alpha_electrons(i)
  unpaired_bitmask(i,2) = unpaired_alpha_electrons(i)
 enddo
 norm_total = 0.d0
 provide mo_two_e_integrals_in_map
 call initialize_density_matrix_osoci
 call bitstring_to_list(inact_bitmask(1,1), occ(1,1), n_occ_beta, N_int)
 print*,''
 print*,''
 print*,'Mulliken spin population analysis'
 accu =0.d0
 do i = 1, nucl_num
  accu += mulliken_spin_densities(i)
  print*,i,nucl_charge(i),mulliken_spin_densities(i)
 enddo
 print*,''
 print*,''
 print*,'DOING FIRST LMCT !!'
 print*,'Threshold_lmct = ',threshold_lmct
 allocate( zero_bitmask(N_int,2) )
  do i = 1, n_inact_orb
   lmct = .True.
   i_hole_osoci = list_inact(i)
   ! First set the current generators to the one of restart
   call check_symetry(i_hole_osoci,thr,test_sym)
   if(.not.test_sym)cycle
   call set_generators_to_generators_restart
   call set_psi_det_to_generators
   call create_restart_and_1h(i_hole_osoci)
   call set_generators_to_psi_det
   call set_bitmask_particl_as_input(reunion_of_bitmask)
   call set_bitmask_hole_as_input(reunion_of_bitmask)
   print*,''
   print*,'--------------------------'
   print*,''
   print*,'Perturbative test for the i_hole_osoci  ',i_hole_osoci
   print*,''
   call is_a_good_candidate(threshold_lmct,is_ok,e_pt2,exit_loop,is_ok_perturbative)
   if(is_ok)then
    print*,''
    print*,''
    print*,''
    print*,'Peturbative test is passed !'
    allocate(dressing_matrix(N_det_generators,N_det_generators))
    dressing_matrix = 0.d0
     do k = 1, N_det_generators
      do l = 1, N_det_generators
        call i_h_j(psi_det_generators(1,1,k),psi_det_generators(1,1,l),N_int,hkl)
        dressing_matrix(k,l) = hkl
      enddo
     enddo
     hkl = dressing_matrix(1,1)
     do k = 1, N_det_generators
       dressing_matrix(k,k) = dressing_matrix(k,k) - hkl
     enddo
     ! Do all the single excitations on top of the CAS and 1h determinants
     call set_bitmask_particl_as_input(reunion_of_bitmask)
     call set_bitmask_hole_as_input(reunion_of_bitmask)
     print*,''
     print*,''
     print*,'Doing all the single excitations ...'
     print*,''
     print*,''
     call all_single(e_pt2)
     call make_s2_eigenfunction
     threshold_davidson = 1.d-6
     soft_touch threshold_davidson 
     call diagonalize_ci
     print*,''
     print*,''
     print*,'Naked matrix'
     do k = 1, N_det_generators
      write(*,'(100(F12.5,X))')dressing_matrix(k,:)
     enddo
     call provide_matrix_dressing(dressing_matrix,n_det_generators,psi_det_generators)

    
     hkl = dressing_matrix(1,1)
     do k = 1, N_det_generators
       dressing_matrix(k,k) = dressing_matrix(k,k) - hkl
     enddo
     print*,''
     print*,'Variational dressing of the matrix'
     do k = 1, N_det_generators
      write(*,'(100(F12.5,X))')dressing_matrix(k,:)
     enddo
     deallocate(dressing_matrix)
   else 
     if(.not.do_it_perturbative)cycle
     if(.not. is_ok_perturbative)cycle
   endif
   call set_intermediate_normalization_lmct_old(norm_tmp,i_hole_osoci)

   do k = 1, N_states
    norm_total(k) += norm_tmp(k)
   enddo
   call update_density_matrix_osoci
 enddo

 if(.True.)then
  print*,''
  print*,'DOING THEN THE MLCT !!'
  print*,'Threshold_mlct = ',threshold_mlct
   lmct = .False.
   do i = 1, n_virt_orb
    i_particl_osoci = list_virt(i)

    print*,'--------------------------'
    ! First set the current generators to the one of restart
    call check_symetry(i_particl_osoci,thr,test_sym)
    if(.not.test_sym)cycle
    call set_generators_to_generators_restart
    call set_psi_det_to_generators
    print*,''
    print*,'--------------------------'
    print*,''
    print*,'Perturbative test for the i_particl_osoci ',i_particl_osoci
    print*,''
    ! Initialize the bitmask to the restart ones
    call initialize_bitmask_to_restart_ones
    ! Impose that only the hole i_hole_osoci can be done
    call modify_bitmasks_for_particl(i_particl_osoci)
    call print_generators_bitmasks_holes
    ! Impose that only the active part can be reached 
    call set_bitmask_hole_as_input(unpaired_bitmask)
    call create_restart_and_1p(i_particl_osoci)
  ! Update the generators 
    call set_generators_to_psi_det
    call set_bitmask_particl_as_input(reunion_of_bitmask)
    call set_bitmask_hole_as_input(reunion_of_bitmask)
  ! so all the single excitation on the new generators 
    call is_a_good_candidate(threshold_mlct,is_ok,e_pt2,exit_loop,is_ok_perturbative)
    if(is_ok)then
      allocate(dressing_matrix(N_det_generators,N_det_generators))
      dressing_matrix = 0.d0
      do k = 1, N_det_generators
       do l = 1, N_det_generators
         call i_h_j(psi_det_generators(1,1,k),psi_det_generators(1,1,l),N_int,hkl)
         dressing_matrix(k,l) = hkl
       enddo
      enddo
      print*,''
      print*,''
      print*,''
      print*,'Peturbative test is passed !'
      print*,''
      print*,''
      print*,'Doing all the single excitations ...'
      print*,''
      print*,''
      call all_single(e_pt2)
      call make_s2_eigenfunction
      threshold_davidson = 1.d-6
      soft_touch threshold_davidson threshold_davidson

      call diagonalize_ci
     
      hkl = dressing_matrix(1,1)
      do k = 1, N_det_generators
        dressing_matrix(k,k) = dressing_matrix(k,k) - hkl
      enddo
      print*,'Naked matrix'
      do k = 1, N_det_generators
       write(*,'(100(F12.5,X))')dressing_matrix(k,:)
      enddo
      call provide_matrix_dressing(dressing_matrix,n_det_generators,psi_det_generators)
      print*,''
      print*,'Variational dressing of the matrix'
      do k = 1, N_det_generators
       write(*,'(100(F12.5,X))')dressing_matrix(k,:)
      enddo
      deallocate(dressing_matrix)
    else
     if(exit_loop)then
      call set_generators_to_generators_restart
      call set_psi_det_to_generators
      exit
     else 
      if(.not.do_it_perturbative)cycle
      if(.not. is_ok_perturbative)cycle
     endif
    endif
    call set_intermediate_normalization_mlct_old(norm_tmp,i_particl_osoci)
    do k = 1, N_states
     norm_total(k) += norm_tmp(k)
    enddo
    call update_density_matrix_osoci
  enddo
 endif
 call set_generators_to_generators_restart
 call set_psi_det_to_generators

   print*,'norm_total = ',norm_total
   norm_total = norm_generators_restart
   norm_total = 1.d0/norm_total
!  call rescale_density_matrix_osoci(norm_total)
   accu = 0.d0
   do i = 1, mo_num
    accu += one_e_dm_mo_alpha_osoci(i,i) + one_e_dm_mo_beta_osoci(i,i)
   enddo
   print*,'accu = ',accu
end
