subroutine Roothaan_Hall_nhSCF

  BEGIN_DOC
  ! non-hermitian SCF
  ! P. SALVADOR, Int. J. Quantum Chem., 2019
  END_DOC

  implicit none

  integer                       :: i, j
  integer                       :: iteration_SCF, dim_DIIS, index_dim_DIIS

  double precision              :: energy_SCF, energy_SCF_previous, delta_energy_SCF
  double precision              :: max_error_DIIS, max_error_DIIS_alpha, max_error_DIIS_beta
  double precision              :: level_shift_save

  double precision, allocatable :: Fock_matrix_DIIS(:,:,:), error_matrix_DIIS(:,:,:)
  double precision, allocatable :: mo_coef_save(:,:)

  logical, external             :: qp_stop


  PROVIDE ao_md5 mo_occ level_shift
  ! character*(32)   :: ao_md5           MD5 key, specific of the |AO| basis
  ! double precision :: mo_occ(mo_num)   |MO| occupation numbers
  ! double precision :: level_shift      Energy shift on the virtual MOs to improve SCF convergence


  allocate( mo_coef_save(ao_num,mo_num) )
  mo_coef_save = 0.d0

  allocate( Fock_matrix_DIIS (ao_num,ao_num,max_dim_DIIS) &
          , error_matrix_DIIS(ao_num,ao_num,max_dim_DIIS) )
  Fock_matrix_DIIS  = 0.d0
  error_matrix_DIIS = 0.d0

  call write_time(6)
  print*,'Energy of the guess = ', SCF_energy
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '  N ', 'Energy  ', 'Energy diff  ',  'DIIS error  ', 'Level shift   '
  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'


  ! initialize energies and density matrices
  energy_SCF_previous = SCF_energy
  delta_energy_SCF    = 1.d0
  iteration_SCF       = 0
  dim_DIIS            = 0
  max_error_DIIS      = 1.d0


  PROVIDE FQS_SQF_matrix_AO Fock_matrix_tc_ao_tot

  ! -----------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------

  do while( ( (max_error_DIIS>threshold_DIIS_nonzero).or.(dabs(delta_energy_SCF)>thresh_SCF) ) &
      .and. ( iteration_SCF < n_it_SCF_max ) )


    ! increment cycle number
    iteration_SCF = iteration_SCF + 1


    ! ---------------------------------------------------------------------
    !
    ! current size of the DIIS space
    dim_DIIS = min(dim_DIIS+1,max_dim_DIIS)

    if( dabs(delta_energy_SCF) > 1.d6 ) then

      ! store Fock and error matrices at each iteration
      index_dim_DIIS = mod(dim_DIIS-1,max_dim_DIIS) + 1
      Fock_matrix_DIIS (1:ao_num,1:ao_num,index_dim_DIIS) = Fock_matrix_tc_ao_tot(1:ao_num,1:ao_num)
      error_matrix_DIIS(1:ao_num,1:ao_num,index_dim_DIIS) = FQS_SQF_matrix_AO(1:ao_num,1:ao_num)

      ! Compute the extrapolated Fock matrix
      call extrapolate_Fock_matrix( error_matrix_DIIS, Fock_matrix_DIIS                  &
                                  , Fock_matrix_tc_ao_tot, size(Fock_matrix_tc_ao_tot,1) &
                                  , iteration_SCF, dim_DIIS )

      Fock_matrix_tc_ao_alpha = Fock_matrix_tc_ao_tot * 0.5d0
      Fock_matrix_tc_ao_beta  = Fock_matrix_tc_ao_tot * 0.5d0
      TOUCH Fock_matrix_tc_ao_alpha Fock_matrix_tc_ao_beta

    endif
    !
    ! ---------------------------------------------------------------------

    mo_coef(1:ao_num,1:mo_num) = Reig_tcFock_mo(1:ao_num,1:mo_num)
    TOUCH MO_coef

    ! Calculate error vectors
    max_error_DIIS = maxval(Abs(FQS_SQF_Matrix_MO))

    ! SCF energy
    energy_SCF       = SCFtc_energy
    Delta_Energy_SCF = energy_SCF - energy_SCF_previous
    if( Delta_Energy_SCF > 0.d0 ) then
      Fock_matrix_tc_ao_tot(1:ao_num,1:ao_num) = Fock_matrix_DIIS(1:ao_num,1:ao_num,index_dim_DIIS)
      Fock_matrix_tc_ao_alpha = Fock_matrix_tc_ao_tot * 0.5d0
      Fock_matrix_tc_ao_beta  = Fock_matrix_tc_ao_tot * 0.5d0
      TOUCH Fock_matrix_tc_ao_alpha Fock_matrix_tc_ao_beta
    endif

    level_shift_save = level_shift
    mo_coef_save(1:ao_num,1:mo_num) = mo_coef(1:ao_num,1:mo_num)

    ! ---------------------------------------------------------------------
    !
    do while(delta_energy_SCF > 0.d0)

      mo_coef(1:ao_num,1:mo_num) = mo_coef_save(1:ao_num,1:ao_num)
      if(level_shift <= .1d0) then
        level_shift = 1.d0
      else
        level_shift = level_shift * 3.0d0
      endif
      TOUCH mo_coef level_shift

      mo_coef(1:ao_num,1:mo_num) = Reig_tcFock_mo(1:ao_num,1:mo_num)
      TOUCH mo_coef

      delta_energy_SCF = SCFtc_energy - energy_SCF_previous
      energy_SCF       = SCFtc_energy
      if(level_shift-level_shift_save > 40.d0) then
        level_shift = level_shift_save * 4.d0
        SOFT_TOUCH level_shift
        exit
      endif
      dim_DIIS = 0

    enddo
    !
    ! ---------------------------------------------------------------------

    level_shift = level_shift * 0.5d0
    SOFT_TOUCH level_shift
    energy_SCF_previous = energy_SCF

    ! print results at the end of each iteration
    write(6,'(I4, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, F16.10, 1X, I3)')  &
      iteration_SCF, energy_SCF, Delta_energy_SCF, max_error_DIIS, level_shift, dim_DIIS

    if(delta_energy_SCF < 0.d0) then
      call save_mos
    endif
    if(qp_stop()) exit

  enddo

  ! -----------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------


  write(6,'(A4, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '====','================','================','================','================'
  write(6,*)
  call write_double(6, Energy_SCF, 'SCF energy')
  call write_time(6)
 
end subroutine Roothaan_Hall_nhSCF

