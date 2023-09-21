!!*********************************************************************
!subroutine e_deriv_cerf(name,fderiv,open,igrad,npt,rho,rhoo,&
!     Ecsrerf,Ecsrerfdr,Ecsrerfdro)
!! Call from dftfun for the calculation of the energy and the first
!! derivatives
!!*********************************************************************

!implicit none

!!input
!logical, intent(in)          :: fderiv,open
!integer, intent(in)          :: npt
!double precision, intent(in) :: rho(1:npt), rhoo(1:npt)

!!output
!character(len=30), intent(out)  :: name
!integer, intent(out)            :: igrad
!double precision, intent(inout) :: Ecsrerf(1:npt), Ecsrerfdr(1:npt)
!double precision, intent(inout) :: Ecsrerfdro(1:npt)

!!local variables
!logical :: fkernel
!double precision, allocatable :: Ecsrerfdrc(:),kernelrho(:),kernelrhoo(:)

!!allocation
!allocate(Ecsrerfdrc(npt),kernelrho(npt),kernelrhoo(npt))

!fkernel = .false.

!call dftfun_ecerf_deriv2(name,fderiv,open,igrad,npt,fkernel,  &
!    rho,rhoo,Ecsrerf,Ecsrerfdrc,Ecsrerfdro,kernelrho,kernelrhoo)
!Ecsrerfdr(:)=Ecsrerfdr(:)+Ecsrerfdrc(:)

!deallocate(Ecsrerfdrc,kernelrho,kernelrhoo)

!end subroutine e_deriv_cerf

!*********************************************************************
!subroutine e_deriv_ldac(name,fderiv,open,igrad,npt,rho,rhoo,&
!    Ec,Ecdr,Ecdro)
!! LDA correlation energy and potentials
!!*********************************************************************
!implicit none

!!input
!logical, intent(in)          :: fderiv,open
!integer, intent(in)          :: npt
!double precision, intent(in) :: rho(1:npt), rhoo(1:npt)

!!output
!character(len=30), intent(out)  :: name
!integer, intent(out)            :: igrad
!double precision, intent(inout) :: Ec(1:npt), Ecdr(1:npt)
!double precision, intent(inout) :: Ecdro(1:npt)

!!local variables
!integer :: i
!double precision :: rs,z,myec,myecrs,myecz,pi
!double precision, allocatable :: Ecdrc(:)


!!allocation
!allocate(Ecdrc(npt))

!!name
!name='LDA correlation'
!igrad=0

!pi = acos(-1d0)

!do i=1,npt

!!Test on the value of rho
!   if (abs(rho(i)).lt.1d-12) then
!        Ec(i)=0d0
!        Ecdrc(i)=0d0
!        Ecdro(i)=0d0
!     else

!        !Closed shell case
!        if (.not.open) then
!           z=0d0
!        else
!           !Open shell case
!           z=rhoo(i)/rho(i)
!        endif

!        rs = (3d0/(4d0*pi*rho(i)))**(1d0/3d0)
!        Ec(i) = myec(rs,z)*rho(i)
!        Ecdrc(i) =  myec(rs,z) -rs/3d0*(myecrs(rs,z))-z*(myecz(rs,z))
!        Ecdro(i) = myecz(rs,z)
!     endif
!enddo

!Ecdr(:)=Ecdr(:)+Ecdrc(:)

!deallocate(Ecdrc)

!end subroutine e_deriv_ldac

!!*********************************************************************
!! subroutine kernel_cerf(spinex,rho,npt,kernelc)
!subroutine kernel_cerf(spinex,kernelc)
!!Call from rpa_response for the calculation of the correlation kernel
!!in the closed-shell case
!! spinex=spin polarizability=0
!!*********************************************************************

!implicit none

!!input
!integer, intent(in)          :: spinex
!! double precision, intent(in) :: rho(1:npt)
!! double precision, intent(in) :: rho_a,rho_b
!
!!output
!double precision, intent(inout) :: kernelc(1:n_points_final_grid)

!!local variables
!character(len=30):: name
!logical          :: fderiv,open,fkernel
!integer          :: igrad
!double precision, allocatable :: rhoo(:),Ecsrerf(:),Ecsrerfdr(:)
!double precision, allocatable :: Ecsrerfdro(:),kernelrho(:),kernelrhoo(:)


!! allocation
!allocate(Ecsrerf(n_points_final_grid),Ecsrerfdr(n_points_final_grid),&
!        &Ecsrerfdro(n_points_final_grid),kernelrho(n_points_final_grid),kernelrhoo(n_points_final_grid))

!fderiv  = .true.
!open    = .false.
!fkernel = .true.

!call dftfun_ecerf_deriv2(name,fderiv,open,igrad,n_points_final_grid,fkernel,  &
!    Ecsrerf,Ecsrerfdr,Ecsrerfdro,kernelrho,kernelrhoo)

!if (spinex.eq.0) then
!  kernelc=kernelrho
!else
!  kernelc=kernelrhoo
!endif

!!write(6,*) 'kernel_cerf', kernelc

!deallocate(rhoo,Ecsrerf,Ecsrerfdr,Ecsrerfdro,kernelrho,kernelrhoo)
!end subroutine kernel_cerf

!***********************************************************************
!subroutine kernel_ldac(spinex,rho,npt,kernelc)
subroutine kernel_ldac(spinex,kernelc)
!***********************************************************************

implicit none
BEGIN_DOC
! LDA correlation kernel
END_DOC

!input
integer, intent(in)         :: spinex
!double precision, intent(in) :: rho(1:npt)

!output
double precision, intent(inout) :: kernelc(n_points_final_grid)

!local variables
integer :: i
double precision :: rs,z,pi,myecrs,myecrs2,myecz2
double precision :: rho, rho_a, rho_b

!write(6,*) 'Calculation of the ldac kernel'

pi = acos(-1d0)
z=0d0

write(*,*) 'n_points_final_grid=', n_points_final_grid

write(*,*) "spinex=",spinex

do i=1,n_points_final_grid
!     rho_a = one_e_dm_and_grad_alpha_in_r(4,i,N_states)
!     rho_b = one_e_dm_and_grad_beta_in_r(4,i,N_states)
     rho_a = one_e_dm_and_grad_alpha_in_r(4,i,1)
     rho_b = one_e_dm_and_grad_beta_in_r(4,i,1)

     rho = rho_a + rho_b

!     call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)

   if (abs(rho).lt.1d-12) then
      kernelc(i) = 0d0
!        write(*,*) 'kernelc(',i,')=', kernelc(i) 
!        write(*,*) 'myecrs=', myecrs(rs,z) 
!        write(*,*) 'myecrs2(',i,')=', myecrs2(rs,z) 
!        write(*,*) '....................................'
   else
      rs = (3d0/(4d0*pi*rho))**(1d0/3d0)
      if (spinex.eq.0) then
         kernelc(i) = -rs*(2d0*(myecrs(rs,z)) - &
              &rs*(myecrs2(rs,z)))/(9d0*rho)

!        write(*,*) 'kernelc(',i,')=', kernelc(i) 
!        write(*,*) 'myecrs=', myecrs(rs,z) 
!        write(*,*) 'myecrs2(',i,')=', myecrs2(rs,z) 
!        write(*,*) '....................................'
        
      else
         kernelc(i) = (myecz2(rs,z))/rho
      endif
   endif
enddo

end subroutine kernel_ldac

!!*********************************************************************
!subroutine dftfun_ecerf_deriv2(name,fderiv,open,igrad,npt,fkernel,  &
!    Ecsrerf,Ecsrerfdr,Ecsrerfdro,kernelrho, kernelrhoo)
! !*********************************************************************
! !
! ! Short range LDA correlation functional complementary to erf
! ! and its first and second derivatives
! ! Author : Elisa Rebolini
! ! Date : 03-2011
! ! From : S. Paziani, S. Moroni, P. Gori-Giorgi and G. Bachelet
! !        PRB,73,155111,(2006)
! !
! !*********************************************************************

!use common_chirs
! implicit none


! ! input
! logical, intent(in)          :: fderiv,open,fkernel
! double precision :: rho, rho_a, rho_b, rho_o, rho_c
! ! output
! character(len=30), intent(out)  :: name
! integer, intent(out)            :: igrad
! double precision, intent(inout) :: Ecsrerf(1:n_points_final_grid), Ecsrerfdr(1:n_points_final_grid)
! double precision, intent(inout) :: Ecsrerfdro(1:n_points_final_grid), kernelrho(1:n_points_final_grid), kernelrhoo(1:n_points_final_grid)

! !double precision, intent(inout) :: Ecsrerfdr2(1:npt),Ecsrerfdro2(1:npt)


! !Local variables
! double precision :: pi,alpha,mu,u,v,Q,a,b,c,d,g1,g2,g3,g4,g5
! double precision :: Acoul,a1,a2,a3,a4,a5
! double precision :: z,phi2,x,rs,b0,C2,C3,C4,C5,cfun4,cfun5
! double precision :: zp,zm,myg2,D2
! double precision :: cd1,cd2,cd3,cd4,cd5,phi8,D3,myec
! double precision :: g,ec1,ec2,eclr

! double precision :: xrs,Qrs,b0rs,grs,C2rs,C3rs,D2rs,D3rs,myecrs
! double precision :: myg2d,cfun4rs,cfun5rs,C4rs,C5rs
! double precision :: a1rs,a2rs,a3rs,a4rs,a5rs,ec1rs,ec2rs,eclrrs

! double precision :: phi2z,xz,uz,vz,Qz,C2z,C3z,zpz,zmz
! double precision :: phi8z,cfun4z,C4z,cfun5z,C5z,myecz
! double precision :: a1z,a2z,a3z,a4z,a5z,eclrz

! double precision :: xrs2,Qrs2,grs2,C2rs2,C3rs2,D2rs2,D3rs2
! double precision :: myg2d2,cfun4rs2,C4rs2,cfun5rs2,C5rs2,myecrs2
! double precision :: a1rs2,a2rs2,a3rs2,a4rs2,a5rs2
! double precision :: ec1rs2,ec2rs2,eclrrs2,ecsrrs2
! double precision :: ct1,ct2,ct3,ct4,ct5,ct6,ct7,ct8

! double precision :: phi2z2,xz2,uz2,vz2,Qz2,C2z2,C3z2,zpz2,zmz2
! double precision :: phi8z2,cfun4z2,C4z2,cfun5z2,C5z2,myecz2
! double precision :: a1z2,a2z2,a3z2,a4z2,a5z2,eclrz2

! double precision :: Ecsrerfdr2,Ecsrerfdro2

! integer :: i

! !name
! name='SR LDA correlation for erf'

! !Constants
! igrad=0
! pi = acos(-1d0)
! alpha = (4d0/(9d0*pi))**(1d0/3d0)
! mu=chirs
! !write(6,*)'mu=',mu
! a   = 5.84605d0
! c   = 3.91744d0
! d   = 3.44851d0
! b   = d -3d0*pi*alpha/(4d0*log(2d0)-4d0)
! g1  = 0.020711d0
! g2  = 0.0819306d0
! g3  = -0.0127713d0
! g4  = 0.00185898d0
! g5  = 0.752411d0
! cd1 = 0.547d0
! cd2 = -0.388d0
! cd3 = 0.676d0
! cd4 = 0.31d0
! cd5 = -4.95d0
! ct1 = -0.776d0
! ct2 = -0.42447200000000007d0
! ct3 = -0.11609309200000006d0
! ct4 = 0.20226528400000005d0
! ct5 = -29.700000000000003d0
! ct6 = -4.138d0
! ct7 = 0.1443049999999999d0
! ct8 = 0.0961d0
! Acoul = (2d0*log(2d0)-2d0)/(pi**2d0)


! !Loop over the grid
! do i=1,n_points_final_grid

!    rho_a = one_e_dm_and_grad_alpha_in_r(4,i,N_states)
!    rho_b = one_e_dm_and_grad_beta_in_r(4,i,N_states)
!    rho = rho_a + rho_b

!    call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)

!
!    !Test on the value of rho
!    if (abs(rho).lt.1d-12) then
!       Ecsrerf(i)=0d0
!       Ecsrerfdr(i)=0d0
!       Ecsrerfdro(i)=0d0
!       kernelrho(i)=0d0
!       kernelrhoo(i)=0d0
!    else
!       
!       
!       !Calculation of eclr per particle

!       !Closed shell case
!       if (.not.open) then
!          z=0d0
!       else
!          !Open shell case
!          z=rho_o/rho
!       endif

!       rs = (3d0/(4d0*pi*rho))**(1d0/3d0)
!       !write(6,'(a,i4)') "i=",i
!       !write(6,'(a,d16.10)') " rho=",rho(i)
!       !write(6,'(a,d16.10)') " z=",z
!       !write(6,'(a,d20.14)') " rs=", rs

!       phi2 = ((1d0+z)**(2d0/3d0) + (1d0-z)**(2d0/3d0))/2d0   !eq(14)

!       x = mu*sqrt(rs)/phi2

!       u = 1d0 + a*x + b*x**2d0 + c*x**3d0

!       v = (1d0+a*x+d*x**2d0)

!       Q = Acoul * &                  !eq(22)
!            log( (1d0 + a*x + b*x**2d0 + c*x**3d0) / (1d0+a*x+d*x**2d0) )

!       b0 = 0.784949d0*rs

!       g =(1d0 + g1*rs + g2*rs**2 + g3*rs**3 + g4*rs**4)*exp(-g5*rs) / 2d0

!       C2 = - 3d0*(1d0-z**2d0)*(g-0.5d0) / (8d0*rs**3d0)      !eq(30)

!       C3 =  - (1d0-z**2d0)*g / (sqrt(2d0*pi)*rs**3d0)       !eq(30)

!       D2 = exp(-cd1*rs)*(cd2*rs + cd3*rs**2d0)/rs**2d0      !eq(33)

!       phi8 = ((1d0+z)**(8d0/3d0) + (1d0-z)**(8d0/3d0))/2d0   !eq(14)

!       D3 = exp(-cd4*rs)*(cd5*rs + rs**2d0) / rs**3d0        !eq(34)

!       !Test on the value of zeta
!       if(abs(z).eq.1.d0) then

!          cfun4 = myg2(rs) -  phi8/(5d0*alpha**2d0*rs**2d0)        !eq(28)

!          cfun5 = myg2(rs)                                         !eq(29)

!       else
!          zp = (2d0/(1d0+z))**(1d0/3d0)

!          zm = (2d0/(1d0-z))**(1d0/3d0)

!          cfun4 = ((1d0 + z)/2d0)**2d0 * myg2(rs*zp) + &           !eq(28)
!               ((1d0 - z)/2d0)**2d0 * myg2(rs*zm) + &
!               (1d0 - z**2d0)*D2 - phi8/(5d0*alpha**2d0*rs**2d0)

!          cfun5 = ((1d0 + z)/2d0)**2d0 * myg2(rs*zp) + &           !eq(29)
!               ((1d0 - z)/2d0)**2d0 * myg2(rs*zm) + (1d0 - z**2d0)*D3

!       endif

!       C4 = - 9d0*cfun4 / (64d0*rs**3d0)                      !eq(30)

!       C5 = - 9d0*cfun5 / (40d0*sqrt(2d0*pi)*rs**3d0)        !eq(30)

!       a1 = 4d0*b0**6d0*C3 + b0**8d0*C5

!       !call ecPW(rs,z,ec,ecd,ecz,ecdd,eczd)

!       !write(6,'(a,d20.14)') " ec=",ec

!       a2 = 4d0*b0**6d0*C2 + b0**8d0*C4 + 6d0*b0**4d0*myec(rs,z)

!       a3 = b0**8d0*C3

!       a4 = b0**8d0*C2 + 4d0*b0**6d0*myec(rs,z)

!       a5 = b0**8d0*myec(rs,z)

!       ec1 = (phi2**3d0*Q + a1*mu**3d0 + a2*mu**4d0 + a3*mu**5d0 &
!            + a4*mu**6d0 + a5*mu**8d0)

!       ec2 = (1d0 + b0**2d0 * mu**2d0)**4d0

!       eclr = ec1 / ec2                                      !eq(26)

!       !************************************************************************
!       !Calculation of Ecsr
!       Ecsrerf(i) = (myec(rs,z) - eclr)*rho
!       !write(6,'(a,d16.10)') " Ecsr=",Ecsrerf(i)

!       !************************************************************************
!       !Test on the value of mu for the first derivatives
!       if (fderiv) then
!          if (mu.eq.0) then
!             Ecsrerfdr(i) =  myec(rs,z)-eclr -rs/3d0*(myecrs(rs,z))-z*(myecz(rs,z))
!             Ecsrerfdro(i) = myecz(rs,z)
!          else
!       !************************************************************************
!       !Calculation of the first derivative of eclr w.r.t. rs

!          xrs = mu / (2d0*phi2*sqrt(rs))

!          Qrs = Acoul * xrs * &
!               ( (a + 2d0*b*x + 3d0*c*x**2d0) / (1d0 + a*x + b*x**2d0 + c*x**3d0) - &
!               (a + 2d0*d*x) / (1d0+a*x+d*x**2d0))

!          b0rs =  0.784949d0

!          grs = exp(-g5*rs) / 2d0 * &
!               ((g1 + 2d0*g2*rs + 3d0*g3*rs**2 + 4d0*g4*rs**3)-&
!               g5*(1d0 + g1*rs + g2*rs**2 + g3*rs**3 + g4*rs**4))

!          C2rs =  - 3d0*(1d0-z**2d0)/8d0 * &
!               (grs*rs - 3d0*(g-0.5d0))/ rs**4d0

!          C3rs = - (1d0-z**2d0) / sqrt(2d0*pi) * &
!               (grs*rs - 3d0*g) / rs**4d0

!          D3rs = -exp(-cd4*rs)*(2*cd5+ rs+cd4*cd5*rs +cd4* rs**2d0) / rs**3d0

!          D2rs = -exp(-cd1*rs)*(cd2 +cd1*cd2*rs + cd1*cd3*rs**2d0)/rs**2d0

!          !Test on the value of zeta
!          if(abs(z).eq.1.d0) then

!             cfun4rs = myg2d(rs) + 2d0*phi8/(5d0*alpha**2d0*rs**3d0)        !eq(28)

!             cfun5rs = myg2d(rs)

!          else

!             cfun4rs = ((1d0 + z)/2d0)**2d0 * zp*myg2d(rs*zp) + &
!                  ((1d0 - z)/2d0)**2d0 * zm*myg2d(rs*zm) + &
!                  (1d0 - z**2d0)*D2rs + 2d0*phi8/(5d0*alpha**2d0*rs**3d0)

!             cfun5rs = ((1d0 + z)/2d0)**2d0 * zp*myg2d(rs*zp) + &
!                  ((1d0 - z)/2d0)**2d0 *zm*myg2d(zm*rs) + (1d0 - z**2d0)*D3rs

!          endif

!          C4rs = - 9d0*(cfun4rs*rs-3d0*cfun4) / (64d0*rs**4d0)

!          C5rs = - 9d0*(cfun5rs*rs-3d0*cfun5) / &
!               (40d0*sqrt(2d0*pi)*rs**4d0)

!          a1rs =  24d0*b0**5d0*b0rs*C3 + 4d0*b0**6d0*C3rs + &
!               (8d0*b0**7d0*b0rs*C5) + b0**8d0*C5rs

!          !write(6,'(a,d16.10)') " myecrs=",myecrs(rs,z)
!          !write(6,'(a,d16.10)') " ecrs=",ecd

!          a2rs = 24d0*b0**5d0*b0rs*C2 + 4d0*b0**6d0*C2rs + &
!               (8d0*b0**7d0*b0rs*C4)  + b0**8d0*C4rs + &
!               (24d0*b0**3d0*b0rs*myec(rs,z))  + 6d0*b0**4d0*myecrs(rs,z)

!          a3rs = 8d0*b0**7d0*b0rs*C3 + b0**8d0*C3rs

!          a4rs = 8d0*b0**7d0*b0rs*C2 + b0**8d0*C2rs + &
!               (24d0*b0**5d0*b0rs*myec(rs,z)) + 4d0*b0**6d0*myecrs(rs,z)

!          a5rs = 8d0*b0**7d0*b0rs*myec(rs,z) + b0**8d0*myecrs(rs,z)

!          ec1rs = phi2**3d0*Qrs + a1rs*mu**3d0 + a2rs*mu**4d0 + a3rs*mu**5d0 &
!               + a4rs*mu**6d0 + a5rs*mu**8d0

!          ec2rs = 8d0*b0*b0rs*mu**2d0*(1d0 + b0**2d0 * mu**2d0)**3d0

!          eclrrs = ec1rs/ec2 - ec1*ec2rs/ec2**2d0

!          !********************************************************************
!          !Calculation of the first derivative of eclr w.r.t. zeta
!          !Test on the value of zeta
!          if (z.gt.1d0-1d-15) then
!             !vrhoo(i)=(2d0*ecz - Eclr + rs*Eclrdr/(rsdr*3d0))/2d0
!             call error('case rhoc=rhoo not implemented ','dftfun_ecerf_deriv2')
!          elseif(z.lt.-1d0+1d-15) then
!             !vrhoo(i)=(2d0*ecz + Eclr - rs*Eclrdr/(rsdr*3d0))/2d0
!             call error('case rhoc=-rhoo not implemented ','dftfun_ecerf_deriv2')
!          else

!             phi2z = (-1d0/(1d0-z)**(1d0/3d0) + 1d0/(1d0+z)**(1d0/3d0))/3d0

!             xz = -mu*sqrt(rs)*phi2z/phi2**2d0

!             uz = (a + 2d0*b*x + 3d0*c*x**2d0)*xz

!             vz = (a + 2d0*d*x)*xz

!             Qz = (2d0*log(2d0)-2d0)/(pi**2d0) * &
!                  (uz/u - vz/v)

!             C2z = 6d0*z*(g-0.5d0) / (8d0*rs**3d0)

!             C3z = 2d0*z*g / (sqrt(2d0*pi)*rs**3d0)

!             phi8z = 4d0*((1d0+z)**(5d0/3d0) -(1d0-z)**(5d0/3d0))/3d0

!             zpz = -2d0**(1d0/3d0) / (3d0*(1d0+z)**(4d0/3d0))

!             zmz = 2d0**(1d0/3d0) / (3d0*(1d0-z)**(4d0/3d0))

!             cfun4z = ((1d0 + z)/2d0) * myg2(rs*zp) + &
!                  ((1d0 + z)/2d0)**2d0 * rs*zpz*myg2d(rs*zp) - &
!                  ((1d0 - z)/2d0) * myg2(rs*zm) + &
!                  ((1d0 - z)/2d0)**2d0 * rs*zmz*myg2d(rs*zm)  &
!                  -2d0*z*D2 - phi8z/(5d0*alpha**2d0*rs**2d0)

!             cfun5z = ((1d0 + z)/2d0) * myg2(rs*zp) + &
!                  ((1d0 + z)/2d0)**2d0 * rs*zpz*myg2d(rs*zp) - &
!                  ((1d0 - z)/2d0) * myg2(rs*zm) + &
!                  ((1d0 - z)/2d0)**2d0 * rs*zmz*myg2d(rs*zm)  &
!                  -2d0*z*D3

!             C4z = - 9d0*cfun4z / (64d0*rs**3d0)

!             C5z = - 9d0*cfun5z / (40d0*sqrt(2d0*pi)*rs**3d0)

!             a1z = 4d0*b0**6d0*C3z + b0**8d0*C5z

!             a2z = 4d0*b0**6d0*C2z + b0**8d0*C4z + 6d0*b0**4d0*myecz(rs,z)

!             a3z = b0**8d0*C3z

!             a4z = b0**8d0*C2z + 4d0*b0**6d0*myecz(rs,z)

!             a5z = b0**8d0*myecz(rs,z)

!             eclrz = (phi2**3d0*Qz + 3d0*phi2z*phi2**2d0*Q + &
!                  a1z*mu**3d0 + a2z*mu**4d0 + a3z*mu**5d0 &
!                  + a4z*mu**6d0 + a5z*mu**8d0) / ec2

!          endif !end of the test on the value of zeta
!          endif !end of the test on the value of mu

!          !*********************************************************************
!          !Calculation of the first derivative of Ecsr w.r.t rho
!          Ecsrerfdr(i) =  myec(rs,z)-eclr -rs/3d0*(myecrs(rs,z)-eclrrs)-z*(myecz(rs,z)-eclrz)
!          !write(6,'(a,d16.10)') " dEdr=", Ecsrerfdr(i)

!          !*********************************************************************
!          !Calculation of the first derivative of Ecsr w.r.t rhoo
!          Ecsrerfdro(i) = myecz(rs,z) - eclrz
!          !write(6,'(a,d16.10)') " dEdro=",Ecsrerfdro(i)


!          if (fkernel) then
!             if (mu.eq.0d0) then

!                kernelrho(i) = -rs*(2d0*(myecrs(rs,z)) - &
!                     rs*(myecrs2(rs,z)))/(9d0*rho)

!                kernelrhoo(i) = (myecz2(rs,z))/rho
!             else
!             !*********************************************************************
!             !Calculation of the second derivative of eclr w.r.t. rs

!             xrs2 = - mu / (4d0*phi2*rs**(3d0/2d0))

!             Qrs2 = Acoul * &
!                  (xrs2 * Qrs/(Acoul*xrs) + xrs**2 * &
!                  ((2d0*b + 6d0*c*x)/(1d0 + a*x + b*x**2 + c*x**3) - &
!                  ((a + 2d0*b*x + 3d0*c*x**2) / (1d0 + a*x + b*x**2 + c*x**3))**2 &
!                  - 2d0*d/(1d0+a*x+d*x**2) + &
!                  ((a + 2d0*d*x) / (1d0+a*x+d*x**2))**2))

!             grs2 = exp(-g5*rs) / 2d0 * &
!                  (-2d0*g5*(g1 + 2d0*g2*rs + 3d0*g3*rs**2 + 4d0*g4*rs**3) +&
!                  g5**2*(1d0 + g1*rs + g2*rs**2 + g3*rs**3 + g4*rs**4) +&
!                  (2d0*g2 + 6d0*g3*rs + 12d0*g4*rs**2))

!             C2rs2 =  - 3d0*(1d0-z**2d0)/8d0 * &
!                  (grs2/rs**3 - 6d0*grs/rs**4 + 12d0*(g-0.5)/rs**5)

!             C3rs2 = - (1d0-z**2d0) / sqrt(2d0*pi) * &
!                  (grs2/rs**3 - 6d0*grs/rs**4 + 12d0*g/rs**5)

!             D2rs2 = exp(-cd1*rs)* &
!                  (ct1 + ct2*rs + ct3*rs**2 + ct4*rs**3)/rs**3

!             cfun4rs2 = ((1d0 + z)/2d0)**2 * zp**2 * myg2d2(rs*zp) + &
!                  ((1d0 - z)/2d0)**2 * zm**2 * myg2d2(rs*zm) + &
!                  (1d0- z**2)*D2rs2 - &
!                  6d0*phi8/(5d0*alpha**2*rs**4)

!             C4rs2 = -9d0*cfun4rs2/(64d0*rs**3) + &
!                  27d0*cfun4rs/(32d0*rs**4) - &
!                  27d0*cfun4/(16d0*rs**5)

!             D3rs2 = exp(-cd4*rs)* &
!                  (ct5 + ct6*rs + ct7*rs**2 + ct8*rs**3)/rs**4

!             cfun5rs2 = ((1d0 + z)/2d0)**2 * zp**2 * myg2d2(rs*zp) + &
!                  ((1d0 - z)/2d0)**2 * zm**2 * myg2d2(rs*zm) + &
!                  (1d0 - z**2)*D3rs2

!             C5rs2 = -9d0/(40d0*sqrt(2d0*pi)) * &
!                  (cfun5rs2/rs**3 - 6d0*cfun5rs/rs**4 +12d0*cfun5/rs**5)

!             a1rs2 = 4d0*b0**6*C3rs2 + 48d0*b0rs*b0**5*C3rs + &
!                  120d0*b0rs**2*b0**4*C3 + b0**8*C5rs2 + &
!                  16d0*b0rs*b0**7*C5rs + 56d0*b0rs**2*b0**6*C5

!             a2rs2 = 4d0*b0**6*C2rs2 + 48d0*b0rs*b0**5*C2rs + &
!                  120d0*b0rs**2*b0**4*C2 + b0**8*C4rs2 + &
!                  16d0*b0rs*b0**7*C4rs + 56d0*b0rs**2*b0**6*C4 + &
!                  6d0*b0**4*myecrs2(rs,z) + 48d0*b0rs*b0**3*myecrs(rs,z) + &
!                  72d0*b0rs**2*b0**2*myec(rs,z)

!             a3rs2 = b0**8*C3rs2 + 16d0*b0rs*b0**7*C3rs + &
!                  56d0*b0rs**2*b0**6*C3

!             a4rs2 = b0**8*C2rs2 + 16d0*b0rs*b0**7*C2rs + &
!                  56d0*b0rs**2*b0**6*C2 + 4d0*b0**6*myecrs2(rs,z) + &
!                  48*b0rs*b0**5*myecrs(rs,z) + 120d0*b0rs**2*b0**4*myec(rs,z)

!             a5rs2 = b0**8*myecrs2(rs,z) + 16d0*b0rs*b0**7*myecrs(rs,z) + &
!                  56d0*b0rs**2*b0**6*myec(rs,z)

!             ec1rs2 = phi2**3*Qrs2 + a1rs2*mu**3 + a2rs2*mu**4 + &
!                  a3rs2*mu**5 + a4rs2*mu**6 + a5rs2*mu**8

!             ec2rs2 = 8d0*b0rs**2*mu**2 * &
!                  ((1d0 +b0**2*mu**2)**3 + 6d0*b0**2*mu**2*(1d0+b0**2*mu**2)**2)

!             eclrrs2 = Ec1rs2/Ec2 - 2d0*Ec1rs*Ec2rs/Ec2**2 - &
!                  Ec1*Ec2rs2/Ec2**2 + 2d0*Ec1*Ec2rs**2/Ec2**3

!             !*********************************************************************
!             !Calculation of the second derivative w.r.t. zeta

!             phi2z2 = -((1d0+z)**(-4d0/3d0) + (1d0-z)**(-4d0/3d0))/9d0

!             xz2 = -mu*sqrt(rs)*(phi2z2/phi2**2 - 2d0*phi2z**2/phi2**3)

!             uz2 = xz2*(a+2d0*b*x+3d0*c*x**2) + xz**2*(2d0*b+6d0*c*x)

!             vz2 = xz2*(a+2d0*d*x) + xz**2*2d0*d

!             Qz2 = (2d0*log(2d0)-2d0)/(pi**2d0) *&
!                  (uz2/u - (uz/u)**2 - vz2/v + (vz/v)**2)

!             C2z2 = 6d0*(g-0.5d0)/(8d0*rs**3)

!             C3z2 = 2d0*g/(sqrt(2d0*pi)*rs**3)

!             zpz2 = 4d0*2d0**(1d0/3d0)/(9d0*(1d0+z)**(7d0/3d0))

!             zmz2 = 4d0*2d0**(1d0/3d0)/(9d0*(1d0-z)**(7d0/3d0))

!             phi8z2 = 20d0*((1d0-z)**(2d0/3d0) + (1d0+z)**(2d0/3d0))/9d0

!             cfun4z2 = ((1d0+z)/2d0)**2* &
!                  (rs*zpz2*myg2d(rs*zp) + rs**2*zpz**2*myg2d2(rs*zp))+&
!                  (1d0+z)*rs*zpz*myg2d(rs*zp) + &
!                  myg2(rs*zp)/2d0 + &
!                  ((1d0-z)/2d0)**2* &
!                  (rs*zmz2*myg2d(rs*zm) + rs**2*zmz**2*myg2d2(rs*zm))-&
!                  (1d0-z)*rs*zmz*myg2d(rs*zm) + &
!                  myg2(rs*zm)/2d0 - &
!                  2d0*D2 - phi8z2/(5d0*alpha**2*rs**2)

!             C4z2 = -9d0*cfun4z2/(64d0*rs**3)

!             cfun5z2 = ((1d0+z)/2d0)**2* &
!                  (rs*zpz2*myg2d(rs*zp) + rs**2*zpz**2*myg2d2(rs*zp))+&
!                  (1d0+z)*rs*zpz*myg2d(rs*zp) + &
!                  myg2(rs*zp)/2d0 + &
!                  ((1d0-z)/2d0)**2* &
!                  (rs*zmz2*myg2d(rs*zm) + rs**2*zmz**2*myg2d2(rs*zm))-&
!                  (1d0-z)*rs*zmz*myg2d(rs*zm) + &
!                  myg2(rs*zm)/2d0 - &
!                  2d0*D3

!             C5z2 = -9d0*cfun5z2/ (40d0*sqrt(2d0*pi)*rs**3d0)

!             a1z2 = 4d0*b0**6d0*C3z2 + b0**8d0*C5z2

!             a2z2 = 4d0*b0**6d0*C2z2 + b0**8d0*C4z2 + &
!                  6d0*b0**4d0*myecz2(rs,z)

!             a3z2 = b0**8d0*C3z2

!             a4z2 = b0**8d0*C2z2 + 4d0*b0**6d0*myecz2(rs,z)

!             a5z2 = b0**8d0*myecz2(rs,z)

!             eclrz2 = (phi2**3*Qz2 + 6d0*phi2z*phi2**2*Qz + &
!                  6d0*phi2z**2*phi2*Q + 3d0*phi2z2*phi2**2*Q + &
!                  a1z2*mu**3d0 + a2z2*mu**4d0 + a3z2*mu**5d0 &
!                  + a4z2*mu**6d0 + a5z2*mu**8d0)/ec2

!             !*******************************************************************
!             !Calculation of the second derivative of Ecsr w.r.t rho when zeta=0
!             !if (spinex.eq.0) then

!                kernelrho(i) = -rs*(2d0*(myecrs(rs,z)-eclrrs) - &
!                     rs*(myecrs2(rs,z) - eclrrs2))/(9d0*rho)

!                !*******************************************************************
!                !Calculation of the second derivative of Ecsr w.r.t. rhoo
!             !else
!                kernelrhoo(i) = (myecz2(rs,z) - eclrz2)/rho

!             !endif !spinex
!          endif !mu
!          !write(6,'(a,d16.10)') "rho=", rho(i)
!          !write(6,'(a,d16.10)') " dEdr2=", kernelrho(i)
!          !write(6,'(a,d16.10)') " dEdro2=", kernelrhoo(i)
!          endif !fkernel
!       endif !fderiv
!    endif !end of the test on the value of rho(i)
! enddo
!end subroutine dftfun_ecerf_deriv2





!*********************************************************************
!Function g"(0,rs,zeta=1) from eq(32)
double precision function myg2(r)

  implicit none

  double precision, intent(in) :: r

  double precision :: gs1,gs2,gs3,pi,alpha

  pi = acos(-1d0)
  alpha = (4d0/(9d0*pi))**(1d0/3d0)
  gs1 = 0.022655d0
  gs2 = 0.4319d0
  gs3 = 0.04d0

  myg2 = 2d0**(5d0/3d0) * (1d0 - gs1*r) / &
       (5d0 * alpha**2d0 * r**2d0 * (1d0 + gs2*r + gs3*r**2d0))

end function myg2



!*********************************************************************
!First derivative of g"
double precision function myg2d(r)

  implicit none

  double precision, intent(in) :: r

  double precision :: pi,alpha,gs1,gs2,gs3,gt1,gt2,gt3

  pi = acos(-1d0)
  alpha = (4d0/(9d0*pi))**(1d0/3d0)
  gs1 = 0.022655d0
  gs2 = 0.4319d0
  gs3 = 0.04d0
  gt1 = -1.273045d0
  gt2 = -0.140430611d0
  gt3 = 0.0027186d0

  myg2d = 2d0**(5d0/3d0)* &
       (-2d0 + gt1*r + gt2*r**2d0 + gt3*r**3d0) / &
       (5d0*alpha**2d0*r**3d0*(1d0 + gs2*r + gs3*r**2d0)**2d0)

end function myg2d


!*********************************************************************
!Second derivatice of g"
double precision function myg2d2(r)

  implicit none

  double precision, intent(in) :: r

  double precision :: pi,alpha,gs2,gs3,gq1,gq2,gq3,gq4,gq5

  pi = acos(-1d0)
  alpha = (4d0/(9d0*pi))**(1d0/3d0)
  gs2 = 0.4319d0
  gs3 = 0.04d0
  gq1 = 6.86509d0
  gq2 = 2.899743153d0
  gq3 = 0.48748674267270004d0
  gq4 = 0.025737795520000002d0
  gq5 = -0.000434976d0

  myg2d2 = 2d0**(5d0/3d0) * (6d0 + gq1*r + gq2*r**2 + &
       gq3*r**3 + gq4*r**4 + gq5*r**5) / &
       (5d0 * alpha**2d0 * r**4d0 * (1d0 + gs2*r + gs3*r**2d0)**3d0)

end function myg2d2





!*******************************************************************
!Function ecPW(rs,zeta) per particle
double precision function myec(rs,z)

  implicit none

  double precision, intent(in) :: rs,z

  double precision :: pi,aaa,f02,f,ec0,ec1,alphac,myGPW

  pi = acos(-1d0)
  aaa = (1d0-log(2d0))/pi**2
  f02 = 4d0/(9d0*(2d0**(1d0/3d0)-1d0))

  f  = ((1d0+z)**(4d0/3d0) + (1d0-z)**(4d0/3d0) - 2d0) / &
       (2d0**(4d0/3d0) -2d0)

  ec0 = myGPW(rs,aaa,0.21370d0,7.5957d0,3.5876d0,1.6382d0,0.49294d0)

  ec1 = myGPW(rs,aaa/2d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,0.62517d0)

  alphac = -myGPW(rs,0.016887d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0)

  myec = ec0 + alphac*f*(1d0-z**4)/f02 + (ec1-ec0)*f*z**4

end function myec

!*******************************************************************
!First derivative of ecPW w.r.t. rs
double precision function myecrs(rs,z)

  implicit none

  double precision, intent(in) :: rs,z

  double precision :: pi,aaa,f02,f,ec0rs,ec1rs,alphacrs,myGPWrs

  pi = acos(-1d0)
  aaa = (1d0-log(2d0))/pi**2
  f02 = 4d0/(9d0*(2d0**(1d0/3d0)-1d0))

  f = ((1d0+z)**(4d0/3d0) + (1d0-z)**(4d0/3d0) - 2d0) / &
       (2d0**(4d0/3d0) -2d0)

  ec0rs = myGPWrs(rs,aaa,0.21370d0,7.5957d0,3.5876d0,1.6382d0,0.49294d0)

  ec1rs = myGPWrs(rs,aaa/2d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,0.62517d0)

  alphacrs = -myGPWrs(rs,0.016887d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0)

  myecrs = ec0rs + alphacrs*f*(1d0-z**4)/f02 + (ec1rs-ec0rs)*f*z**4

end function myecrs

!*******************************************************************
!Second derivative of ecPW w.r.t rs

double precision function myecrs2(rs,z)

  implicit none

  double precision, intent(in) :: rs,z

  double precision :: pi,aaa,f02,f,ec0rs2,ec1rs2,alphacrs2,myGPWrs2
  double precision :: ec,ecd,ecz,ecdd,eczd

  pi = acos(-1d0)
  aaa = (1d0-log(2d0))/pi**2
  f02 = 4d0/(9d0*(2d0**(1d0/3d0)-1d0))

  f = ((1d0+z)**(4d0/3d0) + (1d0-z)**(4d0/3d0) - 2d0) / &
       (2d0**(4d0/3d0) -2d0)

  ec0rs2 = myGPWrs2(rs,aaa,0.21370d0,7.5957d0,3.5876d0,1.6382d0,0.49294d0)

  ec1rs2 = myGPWrs2(rs,aaa/2d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,0.62517d0)

  alphacrs2 = -myGPWrs2(rs,0.016887d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0)

  myecrs2 = ec0rs2 + alphacrs2*f*(1d0-z**4)/f02 + (ec1rs2-ec0rs2)*f*z**4

  !write(6,'(a,d16.10)') " myecrs2=", myecrs2
  !call ecPW(rs,z,ec,ecd,ecz,ecdd,eczd)
  !write(6,'(a,d16.10)') " ecrs2=", ecdd

end function myecrs2

!*******************************************************************
!First derivative of ecPW w.r.t. z
double precision function myecz(rs,z)

  implicit none

  double precision, intent(in) :: rs,z

  double precision ::  pi,aaa,f02,f,fz,ec0,ec1,alphac,myGPW

  pi = acos(-1d0)
  aaa = (1d0-log(2d0))/pi**2
  f02 = 4d0/(9d0*(2d0**(1d0/3d0)-1d0))

  f = ((1d0+z)**(4d0/3d0) + (1d0-z)**(4d0/3d0) - 2d0) / &
       (2d0**(4d0/3d0) -2d0)

  fz = 4d0*((1d0+z)**(1d0/3d0) - (1d0-z)**(1d0/3d0)) / &
       (3d0*(2d0**(4d0/3d0) -2d0))

  ec0 = myGPW(rs,aaa,0.21370d0,7.5957d0,3.5876d0,1.6382d0,0.49294d0)

  ec1 = myGPW(rs,aaa/2d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,0.62517d0)

  alphac = -myGPW(rs,0.016887d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0)

  myecz = alphac*(fz*(1d0-z**4) - 4d0*f*z**3)/f02 + &
       (ec1-ec0)*(fz*z**4 + 4d0*f*z**3)

end function myecz

!*******************************************************************
!Second derivative of ecPW w.r.t z
double precision function myecz2(rs,z)

implicit none

  double precision, intent(in) :: rs,z

  double precision ::  pi,aaa,f02,f,fz,fz2,ec0,ec1,alphac,myGPW

  pi = acos(-1d0)
  aaa = (1d0-log(2d0))/pi**2
  f02 = 4d0/(9d0*(2d0**(1d0/3d0)-1d0))

  f = ((1d0+z)**(4d0/3d0) + (1d0-z)**(4d0/3d0) - 2d0) / &
       (2d0**(4d0/3d0) -2d0)

  fz = 4d0*((1d0+z)**(1d0/3d0) - (1d0-z)**(1d0/3d0)) / &
       (3d0*(2d0**(4d0/3d0) -2d0))

  fz2 = 4d0*((1d0+z)**(-2d0/3d0) + (1d0-z)**(-2d0/3d0)) / &
       (9d0*(2d0**(4d0/3d0) -2d0))

  ec0 = myGPW(rs,aaa,0.21370d0,7.5957d0,3.5876d0,1.6382d0,0.49294d0)

  ec1 = myGPW(rs,aaa/2d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,0.62517d0)

  alphac = -myGPW(rs,0.016887d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0)

  myecz2 = alphac*(fz2*(1d0-z**4) - 8d0*fz*z**3-12d0*f*z**2)/f02 + &
       (ec1-ec0)*(fz2*z**4 + 8d0*f*z**3 + 12d0*f*z**2)

end function myecz2




!*******************************************************************
!Function GPW(rs,A,alpha1,beta1,beta2,beta3,beta4,p)
double precision function myGPW(rs,A,a1,b1,b2,b3,b4)

  implicit none

  double precision, intent(in) :: rs

  double precision :: A,a1,b1,b2,b3,b4,u,v

  u = b1*rs**(1d0/2d0) + b2*rs + b3*rs**(3d0/2d0) + b4*rs**2

  v = 1d0 + 1d0/(2d0*A*u)

  myGPW = -2d0*A * (1d0 + a1*rs) * log(v)

  !write(6,'(a,e20.14)') " myGPW=", myGPW
end function myGPW

!*******************************************************************
!First derivative of GPW w.r.t. rs
double precision function myGPWrs(rs,A,a1,b1,b2,b3,b4)

  implicit none

  double precision, intent(in) :: rs

  double precision :: A,a1,b1,b2,b3,b4,u,v,urs,vrs

  double precision :: dsqrt_rs

  u = b1*rs**(1d0/2d0) + b2*rs + b3*rs**(3d0/2d0) + b4*rs**2

  dsqrt_rs=sqrt(rs)

  urs = b2 + b1/(2d0*dsqrt_rs) + 3d0*b3*dsqrt_rs/2d0 + 2d0*b4*rs

  v = 1d0 + 1d0 / (2d0*A*u)

  vrs = - urs / (2d0*A*u**2)

  myGPWrs = -2d0*A*(a1*log(v) + (1+a1*rs)*vrs/v)

end function myGPWrs


!*******************************************************************
!Second derivative of GPW w.r.t. rs
double precision function myGPWrs2(rs,A,a1,b1,b2,b3,b4)

  implicit none

  double precision, intent(in) :: rs

  double precision :: A,a1,b1,b2,b3,b4,u,v,urs,vrs,urs2,vrs2
  double precision :: G,Gd,Gdd

  double precision :: dsqrt_rs

  u = b1*rs**(1d0/2d0) + b2*rs + b3*rs**(3d0/2d0) + b4*rs**2

  dsqrt_rs=sqrt(rs)

  urs = b2 + b1/(2d0*dsqrt_rs) + 3d0*b3*dsqrt_rs/2d0 + 2d0*b4*rs

  urs2 = -b1*rs**(-3d0/2d0)/4d0 + 3d0*b3/(4d0*dsqrt_rs) + 2d0*b4

  v = 1d0 + 1d0 / (2d0*A*u)

  vrs = -urs / (2d0*A*u**2)

  vrs2 = -(urs2/u**2 - 2d0*urs**2/u**3)/(2d0*A)

  myGPWrs2 = -2d0*A*(2d0*a1*vrs/v + (1d0+a1*rs)*(vrs2/v - vrs**2/v**2))

  !call GPW(rs,A,a1,b1,b2,b3,b4,G,Gd,Gdd)

  !write(6,'(a,e20.14)') " myGPWrs2=", myGPWrs2
  !write(6,'(a,e20.14)') " GPWrs2=", Gdd
end function myGPWrs2
