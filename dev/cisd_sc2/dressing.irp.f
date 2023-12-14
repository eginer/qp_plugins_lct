subroutine get_dressing(i,a,j,b,sigma,tau,delta)
 implicit none
 BEGIN_DOC
 ! get_dressing(i,a,j,b,sigma,tau) returns the diagonal dressing for a double excitation 
 !
 ! of type a^dagger(a,sigma) a(i,sigma) a^dagger(b,tau) a(j,tau)
 END_DOC
 integer, intent(in) :: i,a,j,b,sigma,tau
 double precision, intent(out) :: delta
 delta = 0.d0
 delta += e_cor_a(a,sigma) + e_cor_a(b,tau)
 delta += e_cor_j(j,tau) + e_cor_j(i,sigma)
 delta += e_cor_ab(a,b,sigma,tau)
 delta += e_cor_ij(i,j,sigma,tau)
 delta += e_cor_jb_direct(b,j,tau)
 delta += e_cor_jb_direct(a,i,sigma)
 delta += e_cor_jb_exchange(a,j,sigma,tau)
 delta += e_cor_jb_exchange(b,i,tau,sigma)
 delta += e_cor_jab(b,a,j,tau,sigma)
 delta += e_cor_jab(a,b,i,sigma,tau)
 delta += e_cor_ijb(b,j,i,tau,sigma)
 delta += e_cor_ijb(a,i,j,sigma,tau)
 delta += ecorr_contrib(a,b,i,j,sigma,tau)

end

subroutine get_dressing_bourrin(det_i, delta)
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 BEGIN_DOC
 ! get_dressing(i,a,j,b,sigma,tau) returns the diagonal dressing for a double excitation 
 !
 ! of type a^dagger(a,sigma) a(i,sigma) a^dagger(b,tau) a(j,tau)
 END_DOC
 integer(bit_kind), intent(in) :: det_i(N_int,2)
 double precision, intent(out) :: delta
 integer(bit_kind), allocatable :: det_double(:,:),det_exc(:,:),det_kc(:,:)
 allocate(det_double(N_int,2),det_exc(N_int,2),det_kc(N_int,2))
 integer :: k,l,c,d,s1,s2, i_ok_1, i_ok_2
 integer :: kk,ll,cc,dd,index_i
 delta = 0.d0
 det_double = det_i
 do s1 = 1, 2
  do kk = 1, n_inact_cisd
   k = list_inact_cisd(kk)
   do cc = 1, n_virt_cisd 
    c = list_virt_cisd(cc)
    det_kc = det_double 
    call do_single_excitation(det_kc,k,c,s1,i_ok_1)
    if(i_ok_1 == -1)cycle 
    do s2 = 1, 2
     do ll = 1, n_inact_cisd
       l = list_inact_cisd(ll)
       do dd = 1, n_virt_cisd 
        d = list_virt_cisd(dd)
        det_exc = det_kc
        call do_single_excitation(det_exc,l,d,s2,i_ok_2)
        if(i_ok_2 == -1)cycle 
        delta += ecorr_contrib(dd,cc,ll,kk,s2,s1)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
end

BEGIN_PROVIDER [ double precision, diag_sc2_dressing, (N_det)]
 implicit none
 integer :: ii,jj,aa,bb,s1,s2,a,b,i,j,idet,degree
 double precision :: delta_bourrin, phase
 integer                        :: exc(0:2,2,2)
 diag_sc2_dressing = 0.d0
 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP  PRIVATE(idet,delta_bourrin) &
 !$OMP  SHARED(N_det, psi_det,diag_sc2_dressing)
 do idet = 2, N_det
  call get_dressing_bourrin(psi_det(1,1,idet),delta_bourrin)
  diag_sc2_dressing(idet) = delta_bourrin
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER 
