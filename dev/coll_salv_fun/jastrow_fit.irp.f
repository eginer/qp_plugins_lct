double precision function j_stg_taylor(x,gam,n)
 implicit none
 BEGIN_DOC
! Taylor expansion of exp(-1./gam * exp(-gam*x))
 END_DOC
 double precision, intent(in)  :: x,gam
 integer, intent(in) :: n
 integer :: i
 j_stg_taylor = 1.d0
 do i = 1, n
  j_stg_taylor += fact_inv(i) * (-1.d0)**i * dexp(-dble(i)*gam*x)
 enddo
end

double precision function int_j_stg_taylor(n,D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
  BEGIN_DOC
  ! Approximates the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-1./delta * exp(-delta (r - D)) ) * (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 integer, intent(in) :: n ! order of the Taylor expansion
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

end
