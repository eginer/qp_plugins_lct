subroutine non_hermit_dav(sze,nstates,leigv,reigv,eigval, hcalc_l, hcalc_r)
 implicit none
 BEGIN_DOC
! you enter with a set of "nstates" BI-ORTHONORMAL guess of size "sze" 
!
! reigv(j,i) = <j|r_i>, leigv(j,i) = <j|l_i> with 
!
! <l_j|r_i> = 0 if i =/ j 
!
! and two external subroutine "hcalc_r" and "hcalc_l" which compute H |w> and H^\dagger |w> 
!
! you get out with "nstates" BI-ORTHONORMAL written in leigv,reigv
!
! and the corresponding REAL eigenvalues eigval(i)
 END_DOC
 integer, intent(in) :: sze, nstates
 double precision, intent(inout) :: leigv(sze,nstates), reigv(sze,nstates)
 double precision, intent(inout) :: eigval(nstates)
 external hcalc_r 
 external hcalc_l
 double precision :: htmp,resr,resl,u_dot_v,u_dot_u
 double precision, allocatable :: leigv_tmp(:,:),reigv_tmp(:,:),leigv_full(:,:),reigv_full(:,:)
 double precision, allocatable :: wkr(:),wkr_schmidt(:),wkl(:),wkl_schmidt(:),hr(:),hl(:)
 double precision, allocatable :: leigv_tmp_big(:,:,:),reigv_tmp_big(:,:,:)
 double precision, allocatable :: h_mat(:,:),h_mat_full(:,:),ei(:),s_mat_full(:,:)
 logical :: check_conv(nstates,2)
 allocate(wkr(sze),wkr_schmidt(sze),hr(sze))
 allocate(wkl(sze),wkl_schmidt(sze),hl(sze))
 allocate(leigv_tmp_big(sze,n_states_diag,nstates),reigv_tmp_big(sze,n_states_diag,nstates))
 allocate(h_mat_full(n_states_diag,n_states_diag),leigv_full(n_states_diag,n_states_diag),reigv_full(n_states_diag,n_states_diag))
 allocate(s_mat_full(n_states_diag,n_states_diag))
 integer :: istate,i,j,itmp,n_real_eigv,iter
 istate = 1
 check_conv = .False.
 iter = 0
 do while (.not.check_conv(istate,1))
  iter += 1
  print*,'*******'
  print*,'iter = ',iter
  do i = 1, sze
   leigv_tmp_big(i,1,istate) = leigv(i,istate)
   reigv_tmp_big(i,1,istate) = reigv(i,istate)
  enddo
  s_mat_full = 0.d0
  s_mat_full(istate,istate) = u_dot_v(leigv_tmp_big(1,1,istate),reigv_tmp_big(1,1,istate),sze)
  h_mat_full = 0.d0
  h_mat_full(istate,istate) = eigval(istate)
  do itmp = 2, n_states_diag
   allocate(h_mat(itmp,itmp),ei(itmp),leigv_tmp(itmp,itmp),reigv_tmp(itmp,itmp))

   !!!!!!!!!!!!!! RIGHT COMPONENT 
   call hcalc_r(hr,reigv_tmp_big(1,itmp-1,istate),1,sze) ! hr = H R
   do i = 1, sze
    wkr(i) = hr(i) - eigval(istate) * reigv_tmp_big(i,itmp-1,istate) ! wkr = H R - E R : residual R 
   enddo
   do i = 1, sze
    print*,'wkr',i,wkr(i)
   enddo
   resr = u_dot_u(wkr,sze) ! Norm of the residual R
   wkr = wkr / dsqrt(resr)
   print*,'resr = ',resr
   check_conv(istate,1) = dabs(resr).lt.threshold_davidson
 
   call bi_ortho_gram_schmidt(reigv_tmp_big(1,1,istate),leigv_tmp_big(1,1,istate),sze,itmp-1,wkr,wkr_schmidt)
   do i = 1, sze
    print*,'wkr_schmidt',i,wkr_schmidt(i)
   enddo
   reigv_tmp_big(1:sze,itmp,istate) = wkr_schmidt(1:sze)

 
   !!!!!!!!!!!!!! LEFT COMPONENT 
   call hcalc_l(hl,leigv_tmp_big(1,itmp-1,istate),1,sze) ! hl = H^\dagger L
   do i = 1, sze
    wkl(i) = hl(i) - eigval(istate) * leigv_tmp_big(i,itmp-1,istate) ! wkl = H^\dagger L - E L : residual L 
   enddo
   do i = 1, sze
    print*,'wkl',i,wkl(i)
   enddo
   resl = u_dot_u(wkl,sze) ! Norm of the residual L
   wkl = wkl / dsqrt(resl)
   print*,'resl = ',resl
   check_conv(istate,2) = dabs(resl).lt.threshold_davidson
 
   call bi_ortho_gram_schmidt(leigv_tmp_big(1,1,istate),reigv_tmp_big(1,1,istate),sze,itmp-1,wkl,wkl_schmidt)
   leigv_tmp_big(1:sze,itmp,istate) = wkl_schmidt(1:sze)

   do i = 1, sze
    print*,'wkl_schmidt',i,wkl_schmidt(i)
   enddo

   do i = 1, itmp
    s_mat_full(i,itmp) = u_dot_v(wkl_schmidt,reigv_tmp_big(1,i,istate),sze)
   enddo
   do i = 1, itmp
    s_mat_full(itmp,i) = u_dot_v(wkr_schmidt,leigv_tmp_big(1,i,istate),sze)
   enddo
   print*,''
   print*,''
   print*,'OVERLAP MATRIX '
   do i = 1, itmp
    write(*,'(100(F16.10,X))')s_mat_full(i,1:itmp)
   enddo
   print*,''
   print*,''
   print*,''
   print*,'Building H mat for Right '
   ! building the new matrix elements !
   call hcalc_r(hr,reigv_tmp_big(1,itmp,istate),1,sze) ! hr = H wkr_schmidt
   do i = 1, itmp
    htmp = u_dot_v(leigv_tmp_big(1,i,istate),hr,sze) ! <l_i| H | wkr_schmidt >
    print*,'i,itmp = ',i,itmp
    print*,'htmp',htmp
    h_mat_full(i,itmp) = htmp
   enddo
   print*,'mat tmp'
   do i = 1, itmp
    write(*,'(100(F10.5,X))')h_mat_full(i,1:itmp)
   enddo

   print*,'Building H mat for Left'
   ! building the new matrix elements !
   call hcalc_l(hl,leigv_tmp_big(1,itmp,istate),1,sze) ! hl = H^\dagger wkl_schmidt
   do i = 1, itmp
    htmp = u_dot_v(reigv_tmp_big(1,i,istate),hl,sze) ! <r_i| H^\dagger | wkl_schmidt >
    print*,'i,itmp = ',i,itmp
    print*,'htmp',htmp
    h_mat_full(itmp,i) = htmp
   enddo
   print*,'mat tmp'
   do i = 1, itmp
    write(*,'(100(F10.5,X))')h_mat_full(i,1:itmp)
   enddo

   !!! DIAGONALIZATION OF THE H MATRIX 
   print*,''
   print*,'H matrix '
   print*,''
   do i = 1, itmp
    do j = 1, itmp
     h_mat(j,i) = h_mat_full(j,i)
    enddo
   enddo
   do i = 1, itmp
    write(*,'(100(F10.5,X))')h_mat(i,:)
   enddo
   call non_hrmt_real_diag(itmp,h_mat,reigv_tmp,leigv_tmp,n_real_eigv,ei)
   do i = 1, itmp
    do j = 1, itmp
     leigv_full(j,i) = leigv_tmp(j,i)
     reigv_full(j,i) = reigv_tmp(j,i)
    enddo
   enddo
   print*,'ei(1) = ',ei(1)
   eigval(istate) = ei(1)
   deallocate(h_mat,ei,reigv_tmp,leigv_tmp)
   if(check_conv(istate,1))then
    print*,'converged !!'
    exit
   endif
  enddo
  do itmp = 1, n_states_diag
   do i = 1, sze
    reigv(i,istate) += reigv_full(itmp,istate) * reigv_tmp_big(i,itmp,istate)
    leigv(i,istate) += leigv_full(itmp,istate) * leigv_tmp_big(i,itmp,istate)
   enddo
  enddo
 enddo

end


