program pouet
 implicit none
 integer :: sze,nstates
 double precision, allocatable :: reigv(:,:),leigv(:,:),eigval(:)
 external hcalc_r_tmp
 external hcalc_l_tmp
 sze = n_mat
 nstates = 1
 allocate(reigv(sze,nstates),leigv(sze,nstates),eigval(nstates))
 reigv = 0.d0
 leigv = 0.d0
 reigv(1,1) = 1.d0
 leigv(1,1) = 1.d0
 eigval(1) = h_non_hermit(1,1)
 call non_hermit_dav(sze,nstates,leigv,reigv,eigval,hcalc_l_tmp,hcalc_r_tmp)
 print*,'eigval_ht = ',eigval_ht(1)
 print*,'eigval    = ',eigval(1)
 integer :: i
 do i = 1, sze
  write(*,'(I4,X,2(F16.10,X))')i,reigv(i,1),reigvec_ht(i,1)
 enddo

end
