
BEGIN_PROVIDER [ logical, mo_two_e_integrals_tcdag_int_in_map ]

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
  integer                          :: n_integrals, rc, kk, m, j1, i1, lmax
  double precision                 :: cpu_1, cpu_2, wall_1, wall_2, integral, wall_0
  double precision                 :: map_mb
  integer(ZMQ_PTR)                 :: zmq_to_qp_run_socket, zmq_socket_pull
  integer(map_size_kind)           :: get_mo_tcdag_int_map_size, mo_tcdag_int_map_size
  integer, external                :: add_task_to_taskserver
  integer, external                :: zmq_set_running
  integer(key_kind)                :: key_max
  integer(map_size_kind)           :: sze
  character(len=:),    allocatable :: task
  integer(key_kind),   allocatable :: buffer_i(:)
  real(integral_kind), allocatable :: buffer_value(:)

  PROVIDE mo_tc_sym_two_e_pot_dag_in_map mo_non_hermit_term mo_two_e_integrals_in_map

  print*, 'Providing the mo tcdag_int integrals'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)

  call map_deinit(mo_integrals_tcdag_int_map)
  call two_e_integrals_index(mo_num,mo_num,mo_num,mo_num,key_max)
  sze = key_max
  call map_init(mo_integrals_tcdag_int_map,sze)
  print*,  'mo map initialized : ', sze


  call new_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'mo_integrals_tcdag_int')

  allocate(character(len=mo_num*12) :: task)
  write(fmt,*) '(', mo_num, '(I5,X,I5,''|''))'
  do l=1,mo_num
    write(task,fmt) (i,l, i=1,l)
    if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task)) == -1) then
      stop 'Unable to add task to server'
    endif
  enddo
  deallocate(task)

  if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
    print *,  irp_here, ': Failed in zmq_set_running'
  endif

  PROVIDE nproc
  !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+1)
      i = omp_get_thread_num()
      if (i==0) then
        call mo_two_e_integrals_tcdag_int_in_map_collector(zmq_socket_pull)
      else
        call mo_two_e_integrals_tcdag_int_in_map_slave_inproc(i)
      endif
  !$OMP END PARALLEL

  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'mo_integrals_tcdag_int')


  print*, 'Sorting the map'
  call map_sort(mo_integrals_tcdag_int_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)
  mo_tcdag_int_map_size = get_mo_tcdag_int_map_size()

  print*, 'mo tcdag_int integrals provided:'
  print*, ' Size of mo tcdag_int map :         ', map_mb(mo_integrals_tcdag_int_map) ,'MB'
  print*, ' Number of mo tcdag_int integrals :', mo_tcdag_int_map_size
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

  mo_two_e_integrals_tcdag_int_in_map = .True.

END_PROVIDER


