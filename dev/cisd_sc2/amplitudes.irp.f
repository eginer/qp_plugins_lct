BEGIN_PROVIDER [ integer, n_virt_cisd]
 implicit none
 n_virt_cisd = mo_num - elec_alpha_num
END_PROVIDER 

BEGIN_PROVIDER [ integer, n_inact_cisd]
 implicit none
 n_inact_cisd = elec_alpha_num - n_core_orb
END_PROVIDER 

 BEGIN_PROVIDER [ integer, list_virt_cisd, (n_virt_cisd)]
&BEGIN_PROVIDER [ integer, list_virt_cisd_reverse, (mo_num)]
 implicit none
 integer :: i
 list_virt_cisd_reverse = -1000
 do i = 1, n_virt_cisd
  list_virt_cisd(i) = i+elec_alpha_num
  list_virt_cisd_reverse(i+elec_alpha_num) = i
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ integer, list_inact_cisd, (n_virt_cisd)]
&BEGIN_PROVIDER [ integer, list_inact_cisd_reverse, (mo_num)]
 implicit none
 integer :: i
 list_inact_cisd_reverse = -1000
 do i = 1, n_inact_cisd
  list_inact_cisd(i) = i+n_core_orb
  list_inact_cisd_reverse(i+n_core_orb) = i
 enddo
END_PROVIDER 



 BEGIN_PROVIDER [integer, amplitude_index, (n_virt_cisd, n_virt_cisd, n_inact_cisd, n_inact_cisd,2,2)]
&BEGIN_PROVIDER [double precision, vij_ab, (n_virt_cisd, n_virt_cisd, n_inact_cisd, n_inact_cisd,2,2)]
&BEGIN_PROVIDER [integer, hf_index]
 implicit none
 integer :: idet,degree,i,j,a,b,s1,s2
 integer :: ii,jj,aa,bb
 integer                        :: exc(0:2,2,2)
 double precision               :: phase,hij,tmp
 amplitude_index = -1000
 integer :: icount
 icount = 0
 tmp = 0.d0
 do idet = 1, N_det
  call get_excitation_degree(ref_bitmask,psi_det(1,1,idet),degree,N_int)
  if(degree==0)then
   hf_index = idet
  endif 
  if(degree.ne.2)cycle
  icount += 1
  call get_excitation(ref_bitmask,psi_det(1,1,idet),exc,degree,phase,N_int)
  call decode_exc(exc,degree,i,a,j,b,s1,s2)
  call i_H_j(ref_bitmask , psi_det(1,1,idet),N_int,hij)
  tmp += psi_coef(idet,1) * hij/psi_coef(1,1)
  aa = list_virt_cisd_reverse(a)
  bb = list_virt_cisd_reverse(b)
  ii = list_inact_cisd_reverse(i)
  jj = list_inact_cisd_reverse(j)
  vij_ab(aa,bb,ii,jj,s1,s2) = hij
  amplitude_index(aa,bb,ii,jj,s1,s2) = idet
 enddo
 print*,'icount =',icount
 print*,'N_det  =',N_det
 print*,'tmp    =',tmp
END_PROVIDER 

BEGIN_PROVIDER [double precision, amplitudes, (n_virt_cisd, n_virt_cisd, n_inact_cisd, n_inact_cisd,2,2)]
 implicit none
 integer :: ii,jj,aa,bb
 integer :: i,j,a,b,s1,s2
 integer :: index_i
 double precision :: inv_hf_coef
 inv_hf_coef = 1.d0/psi_coef(hf_index,1)
 amplitudes = 0.d0
 do s2 = 1, 2
  do s1 = 1, 2
   do ii = 1, n_inact_cisd
     do jj = 1, n_inact_cisd
       do aa = 1, n_virt_cisd
         do bb = 1, n_virt_cisd
           index_i = amplitude_index(aa,bb,ii,jj,s1,s2)
           if(index_i.lt.0)cycle
           amplitudes(aa,bb,ii,jj,s1,s2) = psi_coef(index_i,1) * inv_hf_coef 
         enddo
       enddo
     enddo
   enddo
  enddo
 enddo
END_PROVIDER 

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
END_PROVIDER 
