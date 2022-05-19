
 BEGIN_PROVIDER [ double precision, fock_tc_reigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_mo, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, eigval_fock_tc_mo, (mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_fock_tc_eigvec_mo, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
! EIGENVECTORS OF FOCK MATRIX ON THE MO BASIS and their OVERLAP
 END_DOC
 integer :: n_real_tc 
 integer :: i,j,k,l
 print*,'********'
 print*,'Fock_matrix_tc_mo_tot'
 do i = 1, mo_num
  write(*,'(100(F15.10,X))')Fock_matrix_tc_mo_tot(i,:)
 enddo
 print*,'********'
 call non_hrmt_real_diag_new(mo_num,Fock_matrix_tc_mo_tot,&
      fock_tc_leigvec_mo,fock_tc_reigvec_mo,& 
      n_real_tc,eigval_fock_tc_mo)
 print*,'Eigenvalues of Fock_matrix_tc_mo_tot'
 do i = 1, mo_num
  print*,i,eigval_fock_tc_mo(i)
 enddo
 overlap_fock_tc_eigvec_mo = 0.d0
 do i = 1, mo_num
  do k = 1, mo_num
   do l = 1, mo_num
    overlap_fock_tc_eigvec_mo(k,i) +=  fock_tc_leigvec_mo(l,k) * fock_tc_reigvec_mo(l,i) 
   enddo
  enddo
 enddo
 double precision :: accu_d, accu_nd
 accu_d = 0.D0
 accu_nd = 0.D0
 print*,'overlap_fock_tc_eigvec_mo = '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')overlap_fock_tc_eigvec_mo(i,:)
  do k = 1, mo_num
  if(i==k)then
   accu_d += overlap_fock_tc_eigvec_mo(k,i)
  else
   accu_nd += dabs(overlap_fock_tc_eigvec_mo(k,i))
  endif
  enddo 
 enddo
 print*,'accu_d  MO = ',accu_d
 print*,'accu_nd MO = ',accu_nd
 double precision, allocatable :: tmp(:,:),p_lo_mo(:,:),p_ro_mo(:,:)
 allocate(tmp(mo_num,mo_num),p_ro_mo(mo_num, mo_num),p_lo_mo(mo_num, mo_num))
 call get_inverse(fock_tc_leigvec_mo,mo_num,mo_num,p_lo_mo,mo_num) !! Matrix from Left eigvec to usual MOs
 tmp = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    tmp(j,i) += p_lo_mo(j,k) * fock_tc_leigvec_mo(k,i)
   enddo
  enddo 
 enddo
 print*,'identity matrix = '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')tmp(i,:)
 enddo
 call get_inverse(fock_tc_reigvec_mo,mo_num,mo_num,p_ro_mo,mo_num) !! Matrix from Left eigvec to usual MOs
 tmp = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    tmp(j,i) += p_lo_mo(j,k) * fock_tc_leigvec_mo(k,i)
   enddo
  enddo 
 enddo
 print*,'identity matrix = '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')tmp(i,:)
 enddo
 
 tmp = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     tmp(j,i) += eigval_fock_tc_mo(l) * p_lo_mo(k,j) * p_ro_mo(l,i) * overlap_fock_tc_eigvec_mo(k,l)
    enddo
   enddo
  enddo
 enddo
 print*,'rebuilding of Fock matrix'
 accu_nd = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   accu_nd += dabs(tmp(j,i)-Fock_matrix_tc_mo_tot(j,i))
  enddo
  write(*,'(100(F15.10,X))')tmp(i,:)
 enddo
 print*,'accu_nd = ',accu_nd
stop
 
 END_PROVIDER 



 BEGIN_PROVIDER [double precision, s_fock_matrix_reigvec, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, s_fock_matrix_leigvec, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, s_fock_matrix_eigval, (mo_num)]
 implicit none
 integer :: n_real_tc
 call non_hrmt_real_diag_new(mo_num,overlap_fock_tc_eigvec_mo,&
      s_fock_matrix_leigvec,s_fock_matrix_reigvec,& 
      n_real_tc,s_fock_matrix_eigval)
 integer :: i,j,k
 double precision, allocatable :: tmp(:,:)
 allocate(tmp(mo_num, mo_num))
 tmp = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    tmp(j,i) += s_fock_matrix_eigval(k) * s_fock_matrix_reigvec(j,k) * s_fock_matrix_leigvec(i,k)
   enddo
  enddo
 enddo
 print*,'Rebuilding overlap matrix'
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')tmp(i,:)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, fock_tc_reigvec_ao, (ao_num, mo_num)]
&BEGIN_PROVIDER [ double precision, fock_tc_leigvec_ao, (ao_num, mo_num)]
&BEGIN_PROVIDER [ double precision, overlap_fock_tc_eigvec_ao, (mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
! EIGENVECTORS OF FOCK MATRIX ON THE AO BASIS and their OVERLAP
!
! THE OVERLAP SHOULD BE THE SAME AS overlap_fock_tc_eigvec_mo
 END_DOC

 double precision :: tmp(ao_num, mo_num)
 call rotate_mo_coef(mo_coef,fock_tc_reigvec_mo,mo_overlap,fock_tc_reigvec_ao)
 call rotate_mo_coef(mo_coef,fock_tc_leigvec_mo,mo_overlap,fock_tc_leigvec_ao)

 integer :: i,k,q,p
 overlap_fock_tc_eigvec_ao = 0.d0  
 do i = 1, mo_num 
  do k = 1, mo_num
   do p = 1, ao_num
    do q = 1, ao_num
     overlap_fock_tc_eigvec_ao(k,i) += fock_tc_leigvec_ao(q,k) * ao_overlap(q,p) * fock_tc_reigvec_ao(p,i) 
    enddo
   enddo
  enddo
 enddo
 print*,'overlap_fock_tc_eigvec_ao = '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')overlap_fock_tc_eigvec_ao(i,:)
 enddo
 END_PROVIDER

