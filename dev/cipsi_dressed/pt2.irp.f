
subroutine print_tc_pt2
  use selection_types
  implicit none
  BEGIN_DOC
! Selected Full Configuration Interaction with Stochastic selection and PT2.
  END_DOC
  integer                        :: i,j,k,ndet
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  type(pt2_type)                 :: pt2_data, pt2_data_err
  logical, external              :: qp_stop
  logical                        :: print_pt2

  double precision :: rss
  double precision, external :: memory_of_double
  PROVIDE H_apply_buffer_allocated distributed_davidson mo_two_e_integrals_in_map

!  call set_psi_coef_for_delta_to_psi_coef
  N_iter = 1
  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  allocate (zeros(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error

  relative_error=PT2_relative_error

  zeros = 0.d0
  pt2_data % pt2   = -huge(1.e0)
  pt2_data % rpt2  = -huge(1.e0)
  pt2_data % overlap= 0.d0
  pt2_data % variance = huge(1.e0)

  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  print_pt2 = .False.
  call print_CI_dressed(ndet, E_tc,norm,pt2_data,print_pt2)

  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  if (N_det > N_det_max) then
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    N_det = N_det_max
    soft_touch N_det psi_det psi_coef
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    print_pt2 = .False.
    call print_CI_dressed(ndet, E_tc,norm,pt2_data,print_pt2)
  endif

  double precision :: correlation_energy_ratio,E_denom,E_tc,norm

  correlation_energy_ratio = 0.d0

  thresh_it_dav  = 5.d-5
  soft_touch thresh_it_dav

  print_pt2 = .True.
      write(*,'(A)')  '--------------------------------------------------------------------------------'


    to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
    to_select = max(N_states_diag, to_select)

    E_denom = E_tc ! TC Energy of the current wave function 
    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(E_denom, pt2_data, pt2_data_err, relative_error,to_select) ! Stochastic PT2 and selection

    N_iter += 1

!    call copy_H_apply_buffer_to_wf()

!    PROVIDE  psi_coef
!    PROVIDE  psi_det
!    PROVIDE  psi_det_sorted
!
    print *,'******'
    print *,'norm = ',norm
    print *,'******'
    call print_CI_dressed(ndet, E_tc,norm,pt2_data,print_pt2)


end
