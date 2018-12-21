program qp_ao_ints
  use omp_lib
  implicit none
  BEGIN_DOC
! Increments a running calculation to compute AO integral_ijkl_r3
  END_DOC
  integer :: i

  call switch_qp_run_to_master

  zmq_context = f77_zmq_ctx_new ()

  ! Set the state of the ZMQ
  zmq_state = 'ao_integral_ijkl_r3'

  ! Provide everything needed
  double precision :: integral_ijkl_r3, ao_bielec_integral_ijkl_r3
  integral_ijkl_r3 = ao_bielec_integral_ijkl_r3(1,1,1,1)

  character*(64) :: state
  call wait_for_state(zmq_state,state)
  do while (state /= 'Stopped')
    !$OMP PARALLEL DEFAULT(PRIVATE) PRIVATE(i)
    i = omp_get_thread_num()
    call ao_bielec_integrals_ijkl_r3_in_map_slave_tcp(i)
    !$OMP END PARALLEL
    call wait_for_state(zmq_state,state)
  enddo

  print *,  'Done'
end

