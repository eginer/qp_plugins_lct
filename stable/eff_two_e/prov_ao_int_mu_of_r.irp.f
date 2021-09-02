use map_module

!! AO Map
!! ======

BEGIN_PROVIDER [ type(map_type), ao_int_mu_of_r_map ]
  implicit none
  BEGIN_DOC
  ! AO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(ao_num,ao_num,ao_num,ao_num,key_max)
  sze = key_max
  call map_init(ao_int_mu_of_r_map,sze)
  print*,  'AO map initialized : ', sze
END_PROVIDER

double precision function get_ao_two_e_int_mu_of_r(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  logical, external              :: ao_two_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map 
  !DIR$ FORCEINLINE
  if (ao_two_e_integral_zero(i,j,k,l)) then
    tmp = 0.d0
  else
    !DIR$ FORCEINLINE
    call two_e_integrals_index(i,j,k,l,idx)
    !DIR$ FORCEINLINE
    call map_get(map,idx,tmp)
  endif
  result = tmp
end

subroutine get_ao_two_e_int_mu_of_rs(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All i are retrieved for j,k,l fixed.
  ! physicist convention : <ij|kl> 
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)

  integer                        :: i
  integer(key_kind)              :: hash
  logical, external              :: ao_one_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map ao_int_mu_of_r_map

  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  double precision :: get_ao_two_e_int_mu_of_r
  do i=1,sze
    out_val(i) = get_ao_two_e_int_mu_of_r(i,j,k,l,ao_int_mu_of_r_map)
  enddo

end


subroutine get_ao_two_e_int_mu_of_rs_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int

  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: tmp
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero
  PROVIDE ao_two_e_integrals_in_map

  non_zero_int = 0
  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do i=1,sze
    integer, external :: ao_l4
    double precision, external :: ao_two_e_integral
    !DIR$ FORCEINLINE
    if (ao_two_e_integral_zero(i,j,k,l)) then
      cycle
    endif
    call two_e_integrals_index(i,j,k,l,hash)
    call map_get(ao_int_mu_of_r_map, hash,tmp)
    if (dabs(tmp) < ao_integrals_threshold) cycle
    non_zero_int = non_zero_int+1
    out_val_index(non_zero_int) = i
    out_val(non_zero_int) = tmp
  enddo

end


subroutine get_ao_two_e_int_mu_of_rs_non_zero_jl(j,l,thresh,sze_max,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  double precision, intent(in)   :: thresh
  integer, intent(in)            :: j,l, sze,sze_max
  real(integral_kind), intent(out) :: out_val(sze_max)
  integer, intent(out)           :: out_val_index(2,sze_max),non_zero_int

  integer                        :: i,k
  integer(key_kind)              :: hash
  double precision               :: tmp
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero

  PROVIDE ao_two_e_integrals_in_map
  non_zero_int = 0
  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
  do k = 1, sze
   do i = 1, sze
     integer, external :: ao_l4
     double precision, external :: ao_two_e_integral
     !DIR$ FORCEINLINE
     if (ao_two_e_integral_zero(i,j,k,l)) then
       cycle
     endif
     call two_e_integrals_index(i,j,k,l,hash)
     call map_get(ao_int_mu_of_r_map, hash,tmp)
     if (dabs(tmp) < thresh ) cycle
     non_zero_int = non_zero_int+1
     out_val_index(1,non_zero_int) = i
     out_val_index(2,non_zero_int) = k
     out_val(non_zero_int) = tmp
   enddo
  enddo

end


subroutine get_ao_two_e_int_mu_of_rs_non_zero_jl_from_list(j,l,thresh,list,n_list,sze_max,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO two-electron integrals from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  double precision, intent(in)   :: thresh
  integer, intent(in)            :: sze_max
  integer, intent(in)            :: j,l, n_list,list(2,sze_max)
  real(integral_kind), intent(out) :: out_val(sze_max)
  integer, intent(out)           :: out_val_index(2,sze_max),non_zero_int

  integer                        :: i,k
  integer(key_kind)              :: hash
  double precision               :: tmp
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero

  PROVIDE ao_two_e_integrals_in_map
  non_zero_int = 0
  if (ao_one_e_integral_zero(j,l)) then
    out_val = 0.d0
    return
  endif

  non_zero_int = 0
 integer :: kk
  do kk = 1, n_list
   k = list(1,kk)
   i = list(2,kk)
   integer, external :: ao_l4
   double precision, external :: ao_two_e_integral
   !DIR$ FORCEINLINE
   if (ao_two_e_integral_zero(i,j,k,l)) then
     cycle
   endif
   call two_e_integrals_index(i,j,k,l,hash)
   call map_get(ao_int_mu_of_r_map, hash,tmp)
   if (dabs(tmp) < thresh ) cycle
   non_zero_int = non_zero_int+1
   out_val_index(1,non_zero_int) = i
   out_val_index(2,non_zero_int) = k
   out_val(non_zero_int) = tmp
  enddo

end




function get_ao_mu_of_r_map_size()
  implicit none
  integer (map_size_kind) :: get_ao_mu_of_r_map_size
  BEGIN_DOC
  ! Returns the number of elements in the AO map
  END_DOC
  get_ao_mu_of_r_map_size = ao_int_mu_of_r_map % n_elements
end

subroutine clear_ao_mu_of_r_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the AO map
  END_DOC
  call map_deinit(ao_int_mu_of_r_map)
  FREE ao_int_mu_of_r_map
end


subroutine insert_into_ao_int_mu_of_r_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)

  call map_append(ao_int_mu_of_r_map, buffer_i, buffer_values, n_integrals)
end


BEGIN_PROVIDER [ logical, ao_int_mu_of_r_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC

  integer                        :: i,j,k,l
  double precision               :: ao_two_e_integral,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  include 'utils/constants.include.F'

  ! For integrals file
  integer(key_kind),allocatable  :: buffer_i(:)
  integer,parameter              :: size_buffer = 1024*64
  real(integral_kind),allocatable :: buffer_value(:)

  integer                        :: n_integrals, rc
  integer                        :: kk, m, j1, i1, lmax
  character*(64)                 :: fmt

  double precision               :: map_mb
  PROVIDE read_ao_two_e_integrals io_ao_two_e_integrals
!  if (read_ao_two_e_integrals) then
!    print*,'Reading the AO integrals'
!    call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
!    print*, 'AO integrals provided'
!    ao_int_mu_of_r_in_map = .True.
!  else

    print*, 'Providing the AO integrals with DFT corrections '
    call wall_time(wall_0)
    call wall_time(wall_1)
    call cpu_time(cpu_1)

    if (.True.) then
      ! Avoid openMP
      integral = ao_two_e_integral(1,1,1,1)
      double precision, allocatable :: ao_integrals(:,:)
      allocate(ao_integrals(ao_num, ao_num))
      call compute_all_ijkl_for_jl_mu_of_r_int(1,1,ao_integrals)
    endif

    integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
    call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'ao_integrals')

    character(len=:), allocatable :: task
    allocate(character(len=ao_num*12) :: task)
    write(fmt,*) '(', ao_num, '(I5,X,I5,''|''))'
    do l=1,ao_num
      write(task,fmt) (i,l, i=1,l)
      integer, external :: add_task_to_taskserver
      if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task)) == -1) then
        stop 'Unable to add task to server'
      endif
    enddo
    deallocate(task)

    integer, external :: zmq_set_running
    if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
      print *,  irp_here, ': Failed in zmq_set_running'
    endif

    PROVIDE nproc
    !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+1)
        i = omp_get_thread_num()
        if (i==0) then
          call ao_int_mu_of_r_integrals_in_map_collector(zmq_socket_pull)
        else
          call ao_int_mu_of_r_in_map_slave_inproc(i)
        endif
    !$OMP END PARALLEL

    call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'ao_integrals')


    print*, 'Sorting the map'
    call map_sort(ao_int_mu_of_r_map)
    call cpu_time(cpu_2)
    call wall_time(wall_2)
    integer(map_size_kind)         :: get_ao_mu_of_r_map_size, ao_map_size
    ao_map_size = get_ao_mu_of_r_map_size()

    print*, 'AO integrals provided:'
    print*, ' Size of AO map :         ', map_mb(ao_int_mu_of_r_map) ,'MB'
    print*, ' Number of AO integrals :', ao_map_size
    print*, ' cpu  time :',cpu_2 - cpu_1, 's'
    print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

    ao_int_mu_of_r_in_map = .True.

!    if (write_ao_two_e_integrals.and.mpi_master) then
!      call ezfio_set_work_empty(.False.)
!      call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
!      call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
!    endif

!  endif

END_PROVIDER
