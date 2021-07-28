subroutine h_non_hermite(v,u,Hmat,a,N_st,sze) 
 implicit none
 BEGIN_DOC
 ! Template of routine for the application of H
 !
 ! Here, it is done with the Hamiltonian matrix 
 !
 ! on the set of determinants of psi_det 
 !
 ! Computes $v = a * H | u \rangle$ 
 !
 END_DOC
 integer, intent(in)              :: N_st,sze
 double precision, intent(in)     :: u(sze,N_st), Hmat(sze,sze), a
 double precision, intent(inout)  :: v(sze,N_st)
 integer :: i,j,k 
 do k = 1, N_st
  do j = 1, sze
   do i = 1, sze
    v(i,k) += a * u(j,k) * Hmat(i,j)
   enddo
  enddo
 enddo
end


subroutine exp_tau_H(u,v,hmat,tau,et,N_st,sze)
 implicit none
 BEGIN_DOC
! realises v = (1 - tau (H - et)) u
 END_DOC
 integer, intent(in) :: N_st,sze
 double precision, intent(in) :: hmat(sze,sze), u(sze,N_st), tau, et
 double precision, intent(out):: v(sze,N_st)
 double precision :: a
 integer :: i,j
 v = (1.d0 + tau * et) * u 
 a = -1.d0 * tau
 call h_non_hermite(v,u,Hmat,a,N_st,sze)
end

double precision function project_phi0(u,Hmat0,N_st,sze)
 implicit none
 integer, intent(in)              :: N_st,sze
 double precision, intent(in)     :: u(sze,N_st), Hmat0(sze)
 integer :: j
 project_phi0 = 0.d0
 do j = 1, sze
  project_phi0 += u(j,1) * Hmat0(j) 
 enddo
 project_phi0 *= 1.d0 / u(1,1)
end

subroutine project_ground(u,v,hmat,e0,N_st,sze)
 implicit none
 integer, intent(in) :: N_st,sze
 double precision, intent(in)   :: hmat(sze,sze)
 double precision, intent(inout):: u(sze)
 double precision, intent(out)  :: v(sze), e0
 double precision :: et, e_expect, u_dot_v , project_phi0, a , tau, delta_e_min 
 double precision :: wall0, wall1
 double precision, allocatable :: vtmp(:),Hmat0(:),delta_e(:)
 integer :: j,i
 integer, allocatable :: iorder(:)
 allocate(vtmp(sze),Hmat0(sze),delta_e(sze),iorder(sze))
 do i = 1, sze
  iorder(i) = i
  delta_e(i) = dabs(Hmat(i,i) - Hmat(1,1))
 enddo
 call dsort(delta_e,iorder,sze)
 print*,'delta_e',delta_e(1), delta_e(2)
 if(delta_e(1).lt.1.d-10)then
  delta_e_min = delta_e(2)
 else
  delta_e_min = delta_e(1)
 endif
 tau = 0.5d0 * delta_e_min
 print*,'tau = ',tau
 do j = 1, sze
  Hmat0(j) = Hmat(1,j)
 enddo
 j = 0
 e_expect = 0.d0
 et = Hmat(1,1) 
 a = 1.d0
 call wall_time(wall0)
 do while(dabs(e_expect-et).gt.threshold_davidson)
  j += 1
  call exp_tau_H(u,v,hmat,tau,et,N_st,sze)
  call normalize(v(1),sze)
  vtmp = 0.d0
  call h_non_hermite(vtmp,v,Hmat,a,N_st,sze) 
  e_expect = u_dot_v(vtmp,v,sze)
  et = project_phi0(v,Hmat0,N_st,sze)
!  print*,'iteration j ',j
!  print*,'et, e_expect, dabs(e_expect-et)'
!  print*,et, e_expect,dabs(e_expect-et) 
  u = v
 enddo
 e0 = e_expect
 call wall_time(wall1)
 print*,'Convergend in iteration ',j
 print*,'Time to project on the ground state ', wall1 - wall0
end
