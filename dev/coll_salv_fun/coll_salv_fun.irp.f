program coll_salv_fun
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
!  print *, 'Hello world'
!  call test_int_f_tilde
! call test_shank
  call test_int_f_special
end

subroutine test_shank
 implicit none
 integer :: i,n,k
 double precision, allocatable :: array(:)
 double precision :: accu,shank_general,test,pi
 pi = dacos(-1.d0)
 n = 6
 accu = 0.d0
 allocate(array(0:n))
 do i = 0, n
  accu += 4.d0 * (-1.d0)**dble(i)/(2.d0*dble(i)+1.d0)
  array(i) = accu
 enddo
 test = shank_general(array,n,n)
 print*,'shank, error = ',test,dabs(test -pi )
 print*,'accu , error = ',accu,dabs(accu -pi )
 print*,'pi           = ',pi


end
