
subroutine pol_modif_center(A_center, B_center, iorder, A_pol, B_pol) 

  BEGIN_DOC
  ! 
  ! Transform the pol centerd on A:
  !       [ \sum_i ax_i (x-x_A)^i ] [ \sum_j ay_j (y-y_A)^j ] [ \sum_k az_k (z-z_A)^k ] 
  ! to a pol centered on B
  !       [ \sum_i bx_i (x-x_B)^i ] [ \sum_j by_j (y-y_B)^j ] [ \sum_k bz_k (z-z_B)^k ] 
  !
  END_DOC

  ! useful for max_dim
  include 'constants.include.F'

  implicit none

  integer,          intent(in)  :: iorder(3)
  double precision, intent(in)  :: A_center(3), B_center(3)
  double precision, intent(in)  :: A_pol(0:max_dim, 3)
  double precision, intent(out) :: B_pol(0:max_dim, 3)

  integer                       :: i, Lmax

  do i = 1, 3
    Lmax = iorder(i)
    call pol_modif_center_x( A_center(i), B_center(i), Lmax, A_pol(0:Lmax, i), B_pol(0:Lmax, i) ) 
  enddo

  return
end subroutine pol_modif_center



subroutine pol_modif_center_x(A_center, B_center, iorder, A_pol, B_pol) 

  BEGIN_DOC
  ! 
  ! Transform the pol centerd on A:
  !       [ \sum_i ax_i (x-x_A)^i ] 
  ! to a pol centered on B
  !       [ \sum_i bx_i (x-x_B)^i ] 
  !
  ! bx_i = \sum_{j=i}^{iorder} ax_j (x_B - x_A)^(j-i) j! / [ i! (j-i)! ]
  !      = \sum_{j=i}^{iorder} ax_j (x_B - x_A)^(j-i) binom_func(j,i)
  !
  END_DOC

  implicit none

  integer,          intent(in)  :: iorder
  double precision, intent(in)  :: A_center, B_center
  double precision, intent(in)  :: A_pol(0:iorder)
  double precision, intent(out) :: B_pol(0:iorder)

  integer                       :: i, j
  double precision              :: fact_tmp, dx

  double precision              :: binom_func

  dx = B_center - A_center

!  do i = 0, iorder
!    fact_tmp = 0.d0
!    do j = i, iorder
!      tmp      = gamma(dble(j+1)) / gamma(dble(j-i+1)) 
!      fact_tmp = fact_tmp + A_pol(j) * tmp * dx**dble(j-i)
!    enddo
!    B_pol(i) = fact_tmp / gamma(dble(i+1))
!  enddo

  do i = 0, iorder
    fact_tmp = 0.d0
    do j = i, iorder
      fact_tmp += A_pol(j) * binom_func(j, i) * dx**dble(j-i)
    enddo
    B_pol(i) = fact_tmp
  enddo

  return
end subroutine pol_modif_center_x


