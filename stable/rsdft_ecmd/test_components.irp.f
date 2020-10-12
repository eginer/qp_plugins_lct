program test_components
  implicit none
  BEGIN_DOC
  !
  END_DOC

  read_wf=.True.
  touch read_wf
 
  integer                       :: nx, ipoint, istate
  double precision              :: xmax, dx
  double precision              :: Wee
  double precision              :: Wee_lr_mu0_5, Wee_lr_mu1, Wee_lr_mu2, Wee_lr_mu5, Wee_lr_mu10
  double precision              :: mu0_5, mu1, mu2, mu5, mu10
  double precision              :: x, r(3)
  double precision              :: rho, rho_a, rho_b, on_top, rho_2

  double precision              :: pi, dthetamax, dtheta, theta
  double precision              :: r_initial, r1(3), r2(3), r12
  double precision              :: mos_array_r1(mo_num), mos_array_r2(mo_num)
  integer                       :: ntheta, itheta
 
 ! nx = 500
 ! xmax = 2.d0
 ! dx = xmax/dble(nx)
 ! x  = dx
 ! r(:) = nucl_coord_transp(:,1)

 pi = dacos(-1.d0)
 dthetamax = pi
 ntheta = 500
 dtheta = 2.d0*dthetamax / dble(ntheta)
 r_initial = 0.50
 theta = -pi
 ! Travelling electron r1
 r1(1) = r_initial*dsin(theta)
 r1(2) = r_initial*dcos(theta)
 r1(3) = 0.d0
 ! Fixed electron r2
 r2(1) = 0.d0
 r2(2) = r_initial
 r2(3) = 0.d0

  mu0_5  = 0.5d0
  mu1    = 1.d0
  mu2    = 2.d0
  mu5    = 5.d0
  mu10   = 10.d0
  istate = 1
!    write(41,*) '#r1(1)     Wee     Wee_lr_mu0_5     Wee_lr_mu1     Wee_lr_mu2     Wee_lr_mu5     Wee_lr_mu10     rho     on_top   rho_2  mosr1 mos_r2'

    write(44,*)'#theta, rho, rho_2, Wee'
   ! do ipoint=1, nx ! vérifier le nom
   do itheta=1, ntheta
    r12 = dsqrt( (r1(1) - r2(1))**2 + (r1(2) - r2(2))**2  )

    ! Variation principle
    Wee = 1.d0/(r12)
    ! RS-DFT
    Wee_lr_mu0_5 = derf(mu0_5*r12)/(r12)
    Wee_lr_mu1   = derf(mu1*r12)/(r12)
    Wee_lr_mu2   = derf(mu2*r12)/(r12)
    Wee_lr_mu5   = derf(mu5*r12)/(r12)
    Wee_lr_mu10  = derf(mu10*r12)/(r12)

    call give_all_mos_at_r(r1,mos_array_r1)
    call give_all_mos_at_r(r2,mos_array_r2)

     
    call dm_dft_alpha_beta_at_r(r1,rho_a,rho_b)
    rho = rho_a + rho_b

    call give_n2_cas(r1,r2,istate,rho_2)
   
   ! write(44,*) r12, Wee, Wee_lr_mu0_5, Wee_lr_mu1, Wee_lr_mu2, Wee_lr_mu5, Wee_lr_mu10, rho, rho_2
    write(44,*) theta, rho, rho_2, Wee
   
    theta += dtheta
    r1(1) = r_initial*dsin(theta)
    r1(2) = r_initial*dcos(theta)

   enddo
end program 
