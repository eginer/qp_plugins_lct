
BEGIN_PROVIDER [ logical, ao_tc_sym_two_e_pot_dag_in_map ]

  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC

  use f77_zmq
  use map_module

  include 'utils/constants.include.F'

  implicit none

  integer, parameter               :: size_buffer = 1024*64

  character*(64)                   :: fmt
  integer                          :: i, j, k, l
  double precision                 :: cpu_1, cpu_2, wall_1, wall_2, integral, wall_0, map_mb
  integer(map_size_kind)           :: ao_eff_pot_dag_map_size
  character(len=:), allocatable    :: task

  integer(map_size_kind)           :: get_ao_tc_sym_two_e_pot_dag_map_size
  double precision                 :: ao_tc_sym_two_e_pot_dag

  integer(ZMQ_PTR)                 :: zmq_to_qp_run_socket, zmq_socket_pull
  integer, external                :: add_task_to_taskserver
  integer, external                :: zmq_set_running

  PROVIDE nproc

  integral = ao_tc_sym_two_e_pot_dag(1,1,1,1)

  print*, 'Providing the ao_tc_sym_two_e_pot_dag_map integrals'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)

  call new_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'ao_tc_sym_two_e_pot_dag')

  allocate(character(len=ao_num*12) :: task)
  write(fmt,*) '(', ao_num, '(I5,X,I5,''|''))'
  do l = 1, ao_num
    write(task,fmt) (i,l, i = 1, l)
    if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task)) == -1) then
      stop 'Unable to add task to server'
    endif
  enddo
  deallocate(task)

  if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
    print *,  irp_here, ': Failed in zmq_set_running'
  endif

  !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+1)
  i = omp_get_thread_num()
  if (i==0) then
    call ao_tc_sym_two_e_pot_dag_in_map_collector(zmq_socket_pull)
  else
    call ao_tc_sym_two_e_pot_dag_in_map_slave_inproc(i)
  endif
  !$OMP END PARALLEL

  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'ao_tc_sym_two_e_pot_dag')

  print*, 'Sorting the map'
  call map_sort(ao_tc_sym_two_e_pot_dag_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)
  ao_eff_pot_dag_map_size = get_ao_tc_sym_two_e_pot_dag_map_size()

  print*, ' AO eff_pot_dag integrals provided:'
  print*, ' Size of AO eff_pot_dag map :         ', map_mb(ao_tc_sym_two_e_pot_dag_map) ,'MB'
  print*, ' Number of AO eff_pot_dag integrals :', ao_eff_pot_dag_map_size
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

  ao_tc_sym_two_e_pot_dag_in_map = .True.

END_PROVIDER
