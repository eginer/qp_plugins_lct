
 BEGIN_PROVIDER [double precision, ecorr_contrib, (n_virt_cisd, n_virt_cisd, n_inact_cisd, n_inact_cisd,2,2)]
&BEGIN_PROVIDER [double precision, ecorr_tot ]
 implicit none
 integer :: ii,jj,aa,bb
 integer :: i,j,a,b,s1,s2
 integer :: index_i
 double precision :: inv_hf_coef,contrib
 inv_hf_coef = 1.d0/psi_coef(hf_index,1)
 ecorr_tot = 0.d0
 ecorr_contrib = 0.d0
 do s2 = 1, 2
  do s1 = 1, 2
   do jj = 1, n_inact_cisd
     do ii = 1, n_inact_cisd
       do bb = 1, n_virt_cisd
         do aa = 1, n_virt_cisd
           index_i = amplitude_index(aa,bb,ii,jj,s1,s2)
           if(index_i.lt.0)cycle
           contrib =  (psi_coef(index_i,1) * inv_hf_coef * vij_ab(aa,bb,ii,jj,s1,s2))
           ecorr_contrib(aa,bb,ii,jj,s1,s2) = contrib 
           ecorr_tot += contrib 
         enddo
       enddo
     enddo
   enddo
  enddo
 enddo
 print*,'ecorr_tot = ',ecorr_tot
END_PROVIDER 


BEGIN_PROVIDER [ double precision, e_cor_a, (n_virt_cisd,2)]
 implicit none
 integer :: k,l,c,s1,sigma,a
 BEGIN_DOC
!e_cor_a(a,sigma) = - sum_(klc)sum_(sigma_1) e_{k,sigma l,tau}^{a,sigma c,tau}
 END_DOC
 e_cor_a = 0.d0
 do s1 = 1, 2
  do sigma = 1, 2
   do k = 1, n_inact_cisd
    do l = 1, n_inact_cisd
     do a = 1, n_virt_cisd 
      do c = 1, n_virt_cisd
       e_cor_a(a,sigma) += - ecorr_contrib(c,a,k,l,s1,sigma)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, e_cor_j, (n_inact_cisd,2)]
 implicit none
 integer :: k,c,d,tau,s1,j,a
 BEGIN_DOC
!e_cor_j(j,tau) = - sum_(kcd) e_{k,sigma j,tau}^{c,sigma d,tau}
 END_DOC
 e_cor_j = 0.d0
 do tau = 1, 2
  do s1= 1, 2
   do j = 1, n_inact_cisd
    do k = 1, n_inact_cisd
     do c = 1, n_virt_cisd
      do d = 1, n_virt_cisd 
       e_cor_j(j,tau) += - ecorr_contrib(d,c,j,k,tau,s1)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, e_cor_ab, (n_virt_cisd,n_virt_cisd,2,2)]
 implicit none
 integer :: k,l,b,a,tau,sigma
 BEGIN_DOC
!e_cor_ab(a,b,sigma,tau) =  sum_(kl) e_{k,sigma l,tau}^{a,sigma b,tau}
 END_DOC
 e_cor_ab = 0.d0
 do tau = 1, 2
  do sigma = 1, 2
   do l = 1, n_inact_cisd
    do k = 1, n_inact_cisd
     do b = 1, n_virt_cisd
      do a = 1, n_virt_cisd 
       e_cor_ab(a,b,sigma,tau) += ecorr_contrib(a,b,k,l,sigma,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 



BEGIN_PROVIDER [ double precision, e_cor_jb_direct, (n_virt_cisd,n_inact_cisd,2)]
 implicit none
 integer :: k,c,b,j,tau,s1
 BEGIN_DOC
!e_cor_jb(j,tau) =  sum_(kc) e_{k,sigma j,tau}^{c,sigma b,tau}
 END_DOC
 e_cor_jb_direct = 0.d0
 do tau = 1, 2
  do s1= 1, 2
   do j = 1, n_inact_cisd
    do k = 1, n_inact_cisd
     do c = 1, n_virt_cisd
      do b = 1, n_virt_cisd 
       e_cor_jb_direct(b,j,tau) +=  ecorr_contrib(b,c,j,k,tau,s1)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, e_cor_jb_exchange, (n_virt_cisd,n_inact_cisd,2,2)]
 implicit none
 integer :: k,d,a,j,tau,sigma
 BEGIN_DOC
!e_cor_jb(j,tau,sigma) =  sum_(kd) e_{k,sigma j,tau}^{a,sigma d,tau}
 END_DOC
 e_cor_jb_exchange = 0.d0
 do tau = 1, 2
  do sigma= 1, 2
   do k = 1, n_inact_cisd
    do j = 1, n_inact_cisd
     do a = 1, n_virt_cisd
      do d = 1, n_virt_cisd 
       e_cor_jb_exchange(a,j,sigma,tau) +=  ecorr_contrib(d,a,j,k,tau,sigma)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 




BEGIN_PROVIDER [ double precision, e_cor_ij, (n_inact_cisd,n_inact_cisd,2,2)]
 implicit none
 integer :: i,j,c,d,sigma,tau
 BEGIN_DOC
!e_cor_ij(i,j,sigma,tau) =  sum_(c,d) e_{i,sigma j,tau}^{c,sigma d,tau}
 END_DOC
 e_cor_ij = 0.d0
 do sigma= 1, 2
  do tau = 1, 2
   do i = 1, n_inact_cisd
    do j = 1, n_inact_cisd
     do d = 1, n_virt_cisd
      do c = 1, n_virt_cisd 
       e_cor_ij(j,i,tau,sigma) += ecorr_contrib(c,d,j,i,sigma,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, e_cor_jab, (n_virt_cisd,n_virt_cisd,n_inact_cisd,2,2)]
 implicit none
 integer :: j,a,b,k,tau,sigma
 BEGIN_DOC
!e_cor_jab(a,b,j,tau,sigma) = - sum_(k) e_{k,sigma j,tau}^{a,sigma b,tau}
 END_DOC
 e_cor_jab = 0.d0
 do tau = 1, 2
  do sigma= 1, 2
   do k = 1, n_inact_cisd
    do j = 1, n_inact_cisd
     do a = 1, n_virt_cisd
      do b = 1, n_virt_cisd 
       e_cor_jab(b,a,j,sigma,tau) +=  - ecorr_contrib(b,a,j,k,tau,sigma)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [ double precision, e_cor_ijb, (n_virt_cisd,n_inact_cisd,n_inact_cisd,2,2)]
 implicit none
 integer :: i,j,c,b,sigma,tau
 BEGIN_DOC
!e_cor_ijb(b,j,i,sigma,tau) =  sum_(c) e_{i,sigma j,tau}^{c,sigma b,tau}
 END_DOC
 e_cor_ijb = 0.d0
 do sigma= 1, 2
  do tau = 1, 2
   do i = 1, n_inact_cisd
    do j = 1, n_inact_cisd
     do b = 1, n_virt_cisd
      do c = 1, n_virt_cisd 
       e_cor_ijb(b,j,i,tau,sigma) += -ecorr_contrib(c,b,j,i,sigma,tau)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

