program pouet
 implicit none
 read_wf = .True.
 touch read_wf
 !call print_matrix_without_pot_dft
  call print_matrix_with_pot 
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

!do i = 1, ao_num 
!  write(*,"(100(F10.5,x))")(potential_deltarho_ecmd_alpha_ao_0(i,:,1) )
!enddo

!print*,' '
!print*,'*********Delta beta********'

!do i = 1, ao_num 
!  write(*,"(100(F10.5,x))")(potential_deltarho_ecmd_beta_ao_0(i,:,1))
!enddo

!print*,' '
!print*,'*********E_c alpha********'
!do i = 1, ao_num 
!  write(*,"(100(F10.5,x))")(potential_e_c_lda_ecmd_alpha_ao_0(i,:,1) )
!enddo

!print*,' '
!print*,'*********E_c beta********'
!do i = 1, ao_num 
!  write(*,"(100(F10.5,x))")(potential_e_c_lda_ecmd_beta_ao_0(i,:,1) )
!enddo



 print*,' '
 print*,' '
 print*,' '
 print*,'*********REF: e_C alpha ********'

 do i = 1, ao_num 
   write(*,"(100(F10.5,x))")(potential_c_alpha_ao_lda(i,:,1))
 enddo

end
