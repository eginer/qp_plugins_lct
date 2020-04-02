double precision function taylor_gauss(alpha,x,n)
 implicit none
 double precision, intent(in) :: alpha,x
 integer, intent(in) :: n
 integer :: i
 taylor_gauss = 1.d0
 do i = 1,n
  taylor_gauss += (-1.d0*x*x*alpha)**dble(i) * fact_inv(i)
 enddo

end

subroutine taylor_gauss_pol(alpha,n,pol)
 implicit none
 double precision, intent(in) :: alpha
 integer, intent(in) :: n
 double precision, intent(out) :: pol(0:n)
 integer :: i
 pol(0) = 1.d0
 do i = 1,n
  pol(i) = (-alpha)**dble(i)*fact_inv(i)
 enddo

end
