subroutine read_fcidump_1_tc
 implicit none
  double precision :: integral
  integer :: i,j,k,l
  logical :: finished
  print*,'Reading the ONE-body integrals from a TC FCIDUMP'
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
  mo_one_e_integrals = 0.d0 
  do while (.True.)
   read(15,*)integral,i,k,j,l
   if(j.ne.0.and.l.ne.0)cycle
   finished = (j==0).and.(l==0).and.(i==0).and.(k==0)
   if(finished)then
    exit 
   endif
   mo_one_e_integrals(i,k) = integral
   mo_one_e_integrals(k,i) = integral
  enddo
!  do i = 1, mo_num
!   do j = 1, mo_num
!    print*,i,j,mo_one_e_integrals(i,j)
!   enddo
!  enddo
  soft_touch mo_one_e_integrals 
  close(unit=15)
end


subroutine read_fcidump_2_tc(array)
 implicit none
 double precision, intent(out) :: array(mo_num, mo_num, mo_num, mo_num)
 integer :: i,j,k,l
 double precision :: integral
 logical :: finished
  print*,'Reading the TWO-body integrals from a TC FCIDUMP'
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
 
 do while (.True.)
  read(15,*)integral,i,k,j,l
  finished = (j==0).and.(l==0)
  if(finished)then
   exit 
  endif
  array(i,j,k,l) = integral
 enddo
end
