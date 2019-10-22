program pouet
 implicit none
 read_wf = .True.
 touch read_wf
 !call print_matrix_without_pot_dft
  call print_matrix_with_pot 
  call print_mu_of_r
end


subroutine print_matrix_without_pot_dft
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,j
 double precision :: hij 
 double precision, allocatable :: Hamil(:,:)
 allocate(Hamil(N_det,N_det))
 print*,'ref_bitmask_energy = ',ref_bitmask_energy
 print*,'N_det = ',N_det
 do i = 1, N_det
  call debug_det(psi_det(1,1,i),N_int)
  do j = 1, N_det
   call i_H_j(psi_det(1,1,i) , psi_det(1,1,j),N_int,hij)
   Hamil(i,j) = hij 
  enddo
 enddo

 print*,'********withou POT*******'

 do i = 1, N_det
   write(*,"(100(F10.5,x))")Hamil(i,:)
 enddo

end


subroutine print_matrix_with_pot
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: i,j,m,n,degree
 integer :: h2,p2,s2,h1,s1,p1
 integer :: exc(0:2,2,2),n_occ_ab(2)
 double precision :: hij,v_ij,phase 
 integer, allocatable           :: occ(:,:)
 double precision, allocatable :: Hamil(:,:),Hamil0(:,:),Hamilpot(:,:)
 allocate(Hamil(N_det,N_det),Hamil0(N_det,N_det),Hamilpot(N_det,N_det))
 allocate(occ(N_int*bit_kind_size,2))
 print*,'ref_bitmask_energy = ',ref_bitmask_energy
 print*,'N_det = ',N_det
 do i = 1, N_det
  call debug_det(psi_det(1,1,i),N_int)
  do j = 1, N_det
   v_ij=0
   call i_H_j(psi_det(1,1,i) , psi_det(1,1,j),N_int,hij)
   call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)

   if (degree .eq. 1 ) then
    call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)

    
    if (s1 .eq. 1 ) then
    !v_ij = phase*potential_c_alpha_mo_sr_pbe_ueg(h1,p1,1) 
     v_ij = phase*(potential_e_c_lda_ecmd_alpha_mo(h1,p1,1)+potential_deltarho_ecmd_alpha_mo(h1,p1,1)) 
    else
    !v_ij = phase*potential_c_beta_mo_sr_pbe_ueg(h1,p1,1) 
     v_ij = phase*(potential_e_c_lda_ecmd_beta_mo(h1,p1,1)+potential_deltarho_ecmd_beta_mo(h1,p1,1))
    endif

   !v_ij = 0.d0 
   elseif (degree .eq. 0 ) then

    call bitstring_to_list_ab(psi_det(1,1,i), occ, n_occ_ab, N_int)
 
    do m = 1, n_occ_ab(1)
    !v_ij +=  potential_c_alpha_mo_sr_pbe_ueg(occ(m,1),occ(m,1),1)
     v_ij +=  potential_e_c_lda_ecmd_alpha_mo(occ(m,1),occ(m,1),1)+potential_deltarho_ecmd_alpha_mo(occ(m,1),occ(m,1),1)
    enddo

    do m = 1, n_occ_ab(2)
    !v_ij +=  potential_c_beta_mo_sr_pbe_ueg(occ(m,2),occ(m,2),1)
     v_ij +=  potential_e_c_lda_ecmd_beta_mo(occ(m,2),occ(m,2),1)+potential_deltarho_ecmd_beta_mo(occ(m,2),occ(m,2),1)
    enddo
   else
    v_ij=0
   endif

   Hamil0(i,j) = hij
   Hamil(i,j)  = hij+v_ij 
   Hamilpot(i,j)  = v_ij 

  enddo
 enddo


 print*,'**************LDA ecmd***********'
 print*,'*********H without pot DFT********'

 do i = 1, N_det
   write(*,"(100(F10.5,x))")Hamil0(i,:)
 enddo


 print*,' '
 print*,'*********H with pot DFT********'

 do i = 1, N_det
   write(*,"(100(F10.5,x))")Hamil(i,:)
 enddo

 print*,' '
 print*,'*********only pot DFT********'

 do i = 1, N_det
   write(*,"(100(F10.5,x))")Hamilpot(i,:)
 enddo


 print*,' '
 print*,' '

 print*,'*********pot ecmd LDA component ao basis ******'
 print*,'*********Delta alpha***************************'

 do i = 1, ao_num 
   write(*,"(100(F10.5,x))")(potential_deltarho_ecmd_alpha_ao_0(i,:,1) )
 enddo
 print*,' '
 print*,'*********Delta beta********'

 do i = 1, ao_num 
   write(*,"(100(F10.5,x))")(potential_deltarho_ecmd_beta_ao_0(i,:,1))
 enddo

 print*,' '
 print*,'*********E_c alpha********'
 do i = 1, ao_num 
   write(*,"(100(F10.5,x))")(potential_e_c_lda_ecmd_alpha_ao_0(i,:,1) )
 enddo

 print*,' '
 print*,'*********E_c beta********'
 do i = 1, ao_num 
   write(*,"(100(F10.5,x))")(potential_e_c_lda_ecmd_beta_ao_0(i,:,1) )
 enddo



 print*,' '
 print*,' '
 print*,' '
 print*,'*********REF: e_C alpha ********'

 do i = 1, ao_num 
   write(*,"(100(F10.5,x))")(potential_c_alpha_ao_lda(i,:,1))
 enddo

end



subroutine print_mu_of_r               
 implicit none
 integer :: i
 double precision :: r1(3),delta,local_potential,two_bod,mu_airbnb
 double precision :: rhoa,rhob,e_c,vc_a,vc_b
 double precision :: aos_array(ao_num) 

 r1(1) = 0.0d0
 r1(2) = 0.0d0
 r1(3) = -3.0d0
 delta= (3.0d0 - r1(3))/10000 

 integer ::k,j,l
 do i = 1, n_act_orb
  do j = 1, n_act_orb
   do k = 1, n_act_orb 
    do l = 1, n_act_orb
     print*,'<l,k|j,i> = ',l,k,j,i
     print*,two_bod_alpha_beta_mo_physicist(l,k,j,i,1)
    enddo
   enddo
  enddo
 enddo
 stop
 do i = 1,10000 
  r1(3) += delta  
 
 !call f_HF_ab(r1,r1,local_potential,two_bod)
  call f_PSI_ab_routine(r1,r1,local_potential,two_bod)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
    local_potential = 1.d+10
  else 
    local_potential = local_potential /  two_bod
  endif
  mu_airbnb =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0

  call dm_dft_alpha_beta_and_all_aos_at_r(r1,rhoa,rhob,aos_array)
  call ec_lda_sr(mu_airbnb,rhoa,rhob,e_c,vc_a,vc_b)

  write(33,*)r1(3),mu_airbnb,rhoa,e_c,vc_a
  write(44,*)r1(3),mu_airbnb,local_potential*two_bod,two_bod
 !write(44,*)r1(3),mu_airbnb,rhob,e_c,vc_b
 
 enddo

end
