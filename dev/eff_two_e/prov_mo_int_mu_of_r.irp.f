
use map_module

BEGIN_PROVIDER [ logical, mo_two_e_int_mu_of_r_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of integrals of the effective two-e potential 
  END_DOC
  
  integer                        :: i,j,k,l
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
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
  
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map 

  ! TODO FOR READ/WRITE 
  PROVIDE read_mo_int_mu_of_r io_mo_int_mu_of_r
  if (read_mo_int_mu_of_r) then
    print*,'Reading the MO mu of r integrals with the effective potential '
    call map_load_from_disk(trim(ezfio_filename)//'/work/mo_ints_mu_of_r',mo_int_mu_of_r_map)
    print*, 'MO mu of r integrals provided'
    mo_two_e_int_mu_of_r_in_map = .True.
    return
  endif
  
  print*, 'Providing the MO mu of r integrals'
  call wall_time(wall_0)
  call wall_time(wall_1)
  call cpu_time(cpu_1)

  integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
  call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'mo_int_mu_of_r')

  character(len=:), allocatable :: task
  allocate(character(len=mo_num*12) :: task)
  write(fmt,*) '(', mo_num, '(I5,X,I5,''|''))'
  do l=1,mo_num
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
        call mo_two_e_int_mu_of_r_in_map_collector(zmq_socket_pull)
      else
        call mo_two_e_int_mu_of_r_in_map_slave_inproc(i)
      endif
  !$OMP END PARALLEL

  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'mo_int_mu_of_r')


  print*, 'Sorting the map'
  call map_sort(mo_int_mu_of_r_map)
  call cpu_time(cpu_2)
  call wall_time(wall_2)
  integer(map_size_kind)         :: get_mo_mu_of_r_map_size, mo_mu_of_r_map_size
  mo_mu_of_r_map_size = get_mo_mu_of_r_map_size()
  
  print*, 'MO mu of r integrals provided:'
  print*, ' Size of MO mu of r map :         ', map_mb(mo_int_mu_of_r_map) ,'MB'
  print*, ' Number of MO mu of r integrals :', mo_mu_of_r_map_size
  print*, ' cpu  time :',cpu_2 - cpu_1, 's'
  print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'
  
  mo_two_e_int_mu_of_r_in_map = .True.

  ! TODO FOR READ/WRITE 
  if (write_mo_int_mu_of_r) then
    call ezfio_set_work_empty(.False.)
    call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints_mu_of_r',mo_int_mu_of_r_map)
    call ezfio_set_eff_two_e_io_mo_int_mu_of_r("Read")
  endif
  
END_PROVIDER
 

!! MO Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_int_mu_of_r_map ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(mo_num,mo_num,mo_num,mo_num,key_max)
  sze = key_max
  call map_init(mo_int_mu_of_r_map,sze)
  print*,  'MO map initialized : ', sze
END_PROVIDER

 BEGIN_PROVIDER [ integer, mo_int_mu_of_r_cache_min ]
&BEGIN_PROVIDER [ integer, mo_int_mu_of_r_cache_max ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the MOs for which the integrals are in the cache
 END_DOC
 mo_int_mu_of_r_cache_min = max(1,mo_num - 63)
 mo_int_mu_of_r_cache_max = mo_num

END_PROVIDER


double precision function get_mo_two_e_int_mu_of_r(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one MO bi-electronic integral from the MO map
  ! (i,k) = 1, (j,l) = 2
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp
  PROVIDE mo_two_e_int_mu_of_r_in_map mo_int_mu_of_r_cache_min
  !DIR$ FORCEINLINE
!  if (mo_overlap_abs(i,k)*mo_overlap_abs(j,l) < mo_integrals_threshold ) then
!    tmp = 0.d0
!  else
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
!  endif
  result = tmp
end



subroutine get_mo_two_e_ints_mu_of_r(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple MO bi-electronic integral from the MO map .
  ! All i are retrieved for j,k,l fixed.
  ! (j,l) = 1 ; (k,:) = 2
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  
  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh
  PROVIDE mo_two_e_int_mu_of_r_in_map mo_int_mu_of_r_map
  thresh = mo_integrals_threshold
  
!  if (mo_overlap_abs(j,l) < thresh) then
!    out_val = 0.d0
!    return
!  endif
  
  double precision :: get_mo_two_e_int_mu_of_r
  do i=1,sze
    out_val(i) = get_mo_two_e_int_mu_of_r(i,j,k,l,mo_int_mu_of_r_map) 
  enddo
  
end

subroutine get_mo_two_e_ints_mu_of_r_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple MO bi-electronic integral from the MO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  real(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int
  
  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh,tmp
  PROVIDE mo_two_e_int_mu_of_r_in_map
  thresh = mo_integrals_threshold
  
  non_zero_int = 0
!  if (mo_overlap_abs(j,l) < thresh) then
!    out_val = 0.d0
!    return
!  endif
 
  non_zero_int = 0
  do i=1,sze
    integer, external :: mo_l4
    double precision, external :: mo_two_e_int_mu_of_r
    !DIR$ FORCEINLINE
!    if (mo_overlap_abs(i,k)*mo_overlap_abs(j,l) < thresh) then
!      cycle
!    endif
    call two_e_integrals_index(i,j,k,l,hash)
    call map_get(mo_int_mu_of_r_map, hash,tmp)
    if (dabs(tmp) < thresh ) cycle
     non_zero_int = non_zero_int+1
     out_val_index(non_zero_int) = i
     out_val(non_zero_int) = tmp
  enddo
  
end


function get_mo_mu_of_r_map_size()
  implicit none
  integer (map_size_kind) :: get_mo_mu_of_r_map_size
  BEGIN_DOC
  ! Returns the number of elements in the MO map
  END_DOC
  get_mo_mu_of_r_map_size = mo_int_mu_of_r_map % n_elements
end

subroutine clear_mo_mu_of_r_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the MO map
  END_DOC
  call map_deinit(mo_int_mu_of_r_map)
  FREE mo_int_mu_of_r_map
end


subroutine insert_into_mo_int_mu_of_r_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  END_DOC
  
  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  
  call map_append(mo_int_mu_of_r_map, buffer_i, buffer_values, n_integrals)
end

