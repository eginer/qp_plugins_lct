
subroutine get_dressed_ci(u0,energies,res)
 implicit none
 double precision, intent(inout) :: u0(N_det) 
 double precision, intent(inout) :: energies(N_states),res
 double precision, allocatable :: u1(:)
 integer :: j,i
 allocate(u1(N_det))
 print*,'ci_energy_dressed = ',ci_energy_dressed(1)
 u0(:) = ci_eigenvectors_dressed(:,1)
 energies(1) = ci_energy_dressed(1)
!!!! DOES NOT DO THE SAME !!!!
! call diagonalize_CI_dressed
 PROVIDE dressing_column_h 
 do j=1,N_states
   do i=1,N_det
     psi_coef(i,j) = CI_eigenvectors_dressed(i,j)
!     print*,'',i,CI_eigenvectors_dressed(i,j)
   enddo
 enddo
 SOFT_TOUCH psi_coef

 u1 = 0.d0
 call htilde_psi_no_store(psi_det,ci_eigenvectors_dressed,n_det,u1)
 u1(:) -= energies(1) * u0(:)
 res = 0.d0
 do i = 1, N_det
  res += u1(i)*u1(i)
 enddo
 res = dsqrt(res)                                                                                                                                                     
 print*,'Norm of the residual vector ', res

end

subroutine smooth_init_tc
 implicit none
 integer :: j,i,k
 double precision :: res,ebefore,thr,hij,delta_E
 double precision, allocatable :: u1(:),u0(:,:),H_jj(:),dressing_vec(:),H_jj_tmp(:),delta(:)
 double precision, allocatable :: energies(:)
 thr = 1.d-10
 allocate(u1(N_det),u0(N_det,N_states_diag))
 allocate(energies(N_states_diag))
 logical :: converged 
  print*,'Smooth initialization '
  integer :: n_init = 5
  double precision :: dx,accu
  dx = 1.d0/dble(n_init)
  accu = 0.d0
  do k = 1, n_init
   print*,''
   print*,'**************'
   print*,'**************'
   print*,'**************'
   print*,'Iteration j ',k
   dressing_column_h *= accu
   soft_touch dressing_column_h
   call get_dressed_ci(u0,energies,res) 
   accu += dx
   delta_E = dabs(ebefore-energies(1))
   ebefore = energies(1)
   if(k.gt.1)then
    print*,'energies,De = ',energies(1),delta_E
   else
    print*,'energies    = ',energies(1)
   endif

!!!!
!   print*,'ci_energy_dressed = ',ci_energy_dressed(1)
!   u0(:,1) = ci_eigenvectors_dressed(:,1)
!   energies(1) = ci_energy_dressed(1)
!!!!!! DOES NOT DO THE SAME !!!!
!!! call diagonalize_CI_dressed
!   PROVIDE dressing_column_h
!   do j=1,N_states
!     do i=1,N_det
!       psi_coef(i,j) = CI_eigenvectors_dressed(i,j)
!     enddo
!   enddo
!   SOFT_TOUCH psi_coef
!  
!   u1 = 0.d0
!   call htilde_psi_no_store(psi_det,ci_eigenvectors_dressed,n_det,u1)
!   u1(:) -= energies(1) * u0(:,1)
!   res = 0.d0
!   do i = 1, N_det
!    res += u1(i)*u1(i)
!   enddo
!   res = dsqrt(res)                                                                                                                             
!   print*,'Norm of the residual vector ', res

  enddo
  print*,'End of Smooth initialization'
end

subroutine get_e_dressed_scf
 implicit none
 integer :: j,i
 double precision :: res,ebefore,thr,hij,delta_E
 double precision, allocatable :: u1(:),u0(:,:),H_jj(:),dressing_vec(:),H_jj_tmp(:),delta(:)
 double precision, allocatable :: energies(:)
 thr = 1.d-06
 allocate(u1(N_det),u0(N_det,N_states_diag))
 allocate(energies(N_states_diag))
 logical :: converged 

! call smooth_init_tc
 res = 1.d0
 j = 0
 u0 = 0.d0
 u0(1:N_det,1) = psi_coef(1:N_det,1)
 energies = 0.d0
 delta_E = 1.d0
 do while(delta_E .gt. thr)
  j += 1
  print*,''
  print*,'**************'
  print*,'**************'
  print*,'**************'
  print*,'Iteration j ',j
  call get_dressed_ci(u0,energies,res)
  print*,'**************'
  delta_E = dabs(ebefore-energies(1))
  if(j.gt.1)then
   print*,'energies,De = ',energies(1),delta_E
  else
   print*,'energies    = ',energies(1)
  print*,'**************'
  print*,'**************'
  print*,'**************'
  endif
  ebefore = energies(1)
!  call routine_save(u0)
 enddo
 print*,'Converged Self consistent TC energy with N_det = ',N_det
 print*,'Converged TC energy = ',energies(1)



end

 BEGIN_PROVIDER [ double precision, ci_energy_dressed_scf, (n_states) ]
&BEGIN_PROVIDER [ double precision, ci_eigenvectors_dressed_scf, (N_det,n_states) ]
 implicit none
 integer :: j,i
 double precision :: res,ebefore,thr,hij,delta_E
 double precision, allocatable :: u1(:),u0(:,:),H_jj(:),dressing_vec(:),H_jj_tmp(:),delta(:)
 double precision, allocatable :: energies(:)
 integer, save :: icount = 0
 thr = 1.d-10
 allocate(u1(N_det),u0(N_det,N_states_diag))
 allocate(energies(N_states_diag))
 logical :: converged 
 res = 1.d0
 j = 0
 u0 = 0.d0
 energies = 0.d0
 delta_E = 1.d0

 !!!! INITIAL STEP 
! u1 = 0.d0
! if(icount == 0)then
!  call set_dressing_column_h_s(psi_coef)
!  call htilde_psi_no_store(psi_det,psi_coef,n_det,u1)
!  energies(1) = ci_electronic_energy(1)
!  u1(:) -= energies(1) * psi_coef(:,1) 
!  icount += 1
! else
!  call set_dressing_column_h_s(ci_eigenvectors_dressed)
!  call htilde_psi_no_store(psi_det,ci_eigenvectors_dressed,n_det,u1)
!  energies(1) = ci_electronic_energy_dressed(1)
!  u1(:) -= energies(1) *ci_eigenvectors_dressed(:,1) 
! endif
 !!!! 
! res = 0.d0
! do i = 1, N_det
!  res += u1(i)*u1(i)
! enddo
! res = dsqrt(res)                                                                                                                                                     
! print*,'Norm of the intial residual vector ', res
 if(N_det .gt. n_det_max_full)then
  print*,'Smooth initialization '
  integer :: n_init = 5
  double precision :: dx,accu
  dx = 1.d0/dble(n_init)
  accu = 0.d0
  do j = 1, n_init
   print*,''
   print*,'**************'
   print*,'**************'
   print*,'**************'
   print*,'Iteration j ',j
   dressing_column_h *= accu
   soft_touch dressing_column_h
   call diagonalize_CI_dressed
   u0 = ci_eigenvectors_dressed
   print*,'ci_energy_dressed = ',ci_energy_dressed(1)
   energies(1) = ci_energy_dressed(1)
 
   print*,'**************'
   delta_E = dabs(ebefore-energies(1))
   if(j.gt.1)then
    print*,'energies,De = ',energies(1),delta_E
   else
    print*,'energies    = ',energies(1)
   print*,'**************'
   print*,'**************'
   print*,'**************'
   endif
   u1 = 0.d0
   call htilde_psi_no_store(psi_det,ci_eigenvectors_dressed,n_det,u1)
   u1(:) -= energies(1) * u0(:,1)
   res = 0.d0
   do i = 1, N_det
    res += u1(i)*u1(i)
   enddo
   res = dsqrt(res)                                                                                                                                                     
   print*,'Norm of the residual vector ', res
   ebefore = energies(1)
   accu += dx
!!  call routine_save(u0)
  enddo
 endif

 do while(delta_E .gt. thr)
  j += 1
  print*,''
  print*,'**************'
  print*,'**************'
  print*,'**************'
  print*,'Iteration j ',j
  dressing_column_h *= 0.1d0
  call diagonalize_CI_dressed
  u0 = ci_eigenvectors_dressed
  print*,'ci_energy_dressed = ',ci_energy_dressed(1)
  energies(1) = ci_energy_dressed(1)

  print*,'**************'
  delta_E = dabs(ebefore-energies(1))
  if(j.gt.1)then
   print*,'energies,De = ',energies(1),delta_E
  else
   print*,'energies    = ',energies(1)
  print*,'**************'
  print*,'**************'
  print*,'**************'
  endif
  u1 = 0.d0
  call htilde_psi_no_store(psi_det,ci_eigenvectors_dressed,n_det,u1)
  u1(:) -= energies(1) * u0(:,1)
  res = 0.d0
  do i = 1, N_det
   res += u1(i)*u1(i)
  enddo
  res = dsqrt(res)                                                                                                                                                     
  print*,'Norm of the residual vector ', res
  ebefore = energies(1)
!  call routine_save(u0)
 enddo
 ci_energy_dressed_scf = ci_energy_dressed
 ci_eigenvectors_dressed_scf = ci_eigenvectors_dressed
 print*,'Converged Self consistent TC energy with N_det = ',N_det
 print*,'Converged TC energy = ',ci_energy_dressed_scf(1)



END_PROVIDER 



