program nh_scf

  BEGIN_DOC 
  ! non-hermitian scf
  END_DOC

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call create_guess
  call orthonormalize_mos
  call run

end


!__________________________________________________________________________________________________
!
subroutine create_guess

  BEGIN_DOC
  !Create a MO guess if no MOs are present in the EZFIO directory
  END_DOC

  implicit none
  logical :: exists

  PROVIDE ezfio_filename
  call ezfio_has_mo_basis_mo_coef(exists)
  if (.not.exists) then
    mo_label = 'Guess'
    if (mo_guess_type == "HCore") then
      mo_coef = ao_ortho_lowdin_coef
      call restore_symmetry(ao_num,mo_num,mo_coef,size(mo_coef,1),1.d-10)
      TOUCH mo_coef
      call mo_as_eigvectors_of_mo_matrix(mo_one_e_integrals,     &
          size(mo_one_e_integrals,1),                            &
          size(mo_one_e_integrals,2),                            &
          mo_label,1,.false.)
      call restore_symmetry(ao_num,mo_num,mo_coef,size(mo_coef,1),1.d-10)
      SOFT_TOUCH mo_coef
    else if (mo_guess_type == "Huckel") then
      call huckel_guess
    else
      print *,  'Unrecognized MO guess type : '//mo_guess_type
      stop 1
    endif
    SOFT_TOUCH mo_label
  endif

end subroutine create_guess
!__________________________________________________________________________________________________


!__________________________________________________________________________________________________
!
subroutine run

  BEGIN_DOC
  ! run non-hermitian SCF calculation
  END_DOC

  implicit none

  mo_label = 'Orthonormalized'

  !thresh_SCF = 1d-2
  !TOUCH thresh_SCF
  !call Roothaan_Hall_SCF
  !call ezfio_set_hartree_fock_energy(SCF_energy)

  !thresh_SCF = 1d-6
  !TOUCH thresh_SCF

  call Roothaan_Hall_nhSCF
  call ezfio_set_hartree_fock_energy(SCF_energy)

end subroutine run
!__________________________________________________________________________________________________

