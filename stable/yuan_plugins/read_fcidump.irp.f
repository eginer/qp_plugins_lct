program read_fcidump
  mo_two_e_integrals_in_map = .True.
  call routine
end
subroutine routine 
  use bitmasks
  use map_module
 implicit none
 double precision :: integral
 integer :: i,j,k,l,n_integrals
 logical :: finished
 open (unit=15, file="FCIDUMP", status='old',    &
              access='sequential', action='read' )
 character*200 :: tmp
 read(15,*)tmp
 print*,tmp
 read(15,*)tmp
 print*,tmp
 read(15,*)tmp
 print*,tmp
 read(15,*)tmp

 integer :: size_buffer 
 integer(key_kind),allocatable  :: buffer_i(:)
 real(integral_kind),allocatable :: buffer_value(:)
 size_buffer = min(ao_num*ao_num*ao_num,16000000)
 provide mo_integrals_map
 allocate(buffer_i(size_buffer),buffer_value(size_buffer))
 n_integrals = 0
 do while (.True.)
  read(15,*)integral,i,k,j,l
  finished = (j==0).and.(l==0)
  if(finished)then
   exit 
  endif
  n_integrals += 1 
  call mo_two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
  buffer_value(n_integrals) = 1.0d0 * integral
!  print*,integral,i,k,j,l
  if (n_integrals == size_buffer) then
    call map_append(mo_integrals_map, buffer_i, buffer_value, n_integrals)
    n_integrals = 0
  endif
 enddo
 print*,'exit main loop'
 print*,n_integrals
  call map_append(mo_integrals_map, buffer_i, buffer_value, n_integrals)
  n_integrals = 0

  integer*8                      :: get_mo_map_size, mo_map_size
  double precision :: get_two_e_integral
  mo_map_size = get_mo_map_size()
  print*,'mo_map_size = ',mo_map_size
  integral = get_two_e_integral(1,1,1,1,mo_integrals_map)
  print*,'<11|11> = ',integral
  call ezfio_set_work_empty(.False.)
  call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
  call ezfio_set_mo_two_e_ints_io_mo_two_e_integrals('Read')


end
