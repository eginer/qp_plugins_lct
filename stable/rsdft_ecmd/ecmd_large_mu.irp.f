 BEGIN_PROVIDER [double precision,integral_on_top,(N_states)]
 implicit none
 BEGIN_DOC
 ! Correlation energy at large mu given by the expression in J. Toulouse, F. Colonna, and A. Savin, Phys. Rev. A 70, 062505 (2004)  and  P. Gori-Giorgi and A. Savin, Phys. Rev. A 73, 032506 (2006)
 ! Ecmd = 2*sqrt(pi) (1 - sqrt(2))/(3 * mu_erf_dft**3) * int(dr) n^(2)(r)
 END_DOC
 integer :: i,j,k,l,istate,t,icouple,jcouple
 integer, allocatable :: order_loc(:)
 double precision :: get_mo_bielec_integral_ijkl_r3
 double precision, allocatable :: E_cor_coupl(:,:)
 allocate(E_cor_coupl(mo_tot_num*mo_tot_num,2),order_loc(mo_tot_num*mo_tot_num))

 double precision, allocatable :: integrals_ij(:,:)
 allocate(integrals_ij(mo_tot_num,mo_tot_num))
 double precision :: wall0,wall1
 integral_on_top = 0.d0
 provide two_bod_alpha_beta_mo_physicist mo_integrals_ijkl_r3_map
 call wall_time(wall0)
   do j = 1, mo_tot_num ! loop over the second electron 
    do i = 1, mo_tot_num ! loop over the first electron 
     E_cor_coupl(icouple,:) = 0.d0
     call get_mo_bielec_integrals_ijkl_r3_ij(i,j,mo_tot_num,integrals_ij,mo_integrals_ijkl_r3_map)
     do l = 1, mo_tot_num
      do k = 1, mo_tot_num
       do istate = 1, N_states
        integral_on_top(istate) += two_bod_alpha_beta_mo_physicist(k,l,i,j,istate) * integrals_ij(k,l)
       enddo
      enddo
     enddo
    enddo
   enddo
 deallocate(integrals_ij)
 call wall_time(wall1)
 print*,'Time to provide integral_on_top   = ',wall1-wall0
 END_PROVIDER

 BEGIN_PROVIDER [double precision, Energy_c_md_on_top, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density.
  ! Ec_md_on_top = (alpha/mu**3) * int n2(r,r) dr  where alpha = sqrt(2pi)*(-2+sqrt(2)) 
 END_DOC
 implicit none 
 integer :: istate
 double precision :: pi,mu
 mu = mu_erf_dft
 pi = 4.d0 * datan(1.d0)
 Energy_c_md_on_top = ((-2.d0+sqrt(2.d0))*sqrt(2.d0*pi)/(3.d0*(mu**3)))*integral_on_top
 END_PROVIDER

