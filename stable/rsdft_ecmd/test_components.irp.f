program test_components
  implicit none
  BEGIN_DOC
  !
  END_DOC

  integer                       :: nx, ipoint
  double precision              :: xmax, dx
  double precision              :: Wee, Wee_lr_mu1, Wee_lr_mu2, mu
  double precision              :: x
  nx = 500
  xmax = 4.d0
  dx = xmax/dble(nx)
  x  = dx

  mu1 = 1.d0
  mu2 = 2.d0

   do ipoint=1, nx ! vérifier le nom
    Wee    = 0.d0 
    Wee_lr_mu1 = 0.d0
    Wee_lr_mu2 = 0.d0
   
    ! Variation principle
    Wee = 1.d0/(x)

    ! RS-DFT
    Wee_lr_mu1 = derf(mu1*x)/(x)
    Wee_lr_mu2 = derf(mu2*x)/(x)

    ! Basis correction
    !

    x +=dx

    write(41,*) x, Wee, Wee_lr_mu1, Wee_lr_mu2
   enddo
end program 
