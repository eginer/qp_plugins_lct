
 BEGIN_PROVIDER [double precision, two_bod_alpha_beta_mo_contracted, (mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num,N_states)]
 implicit none
 BEGIN_DOC
 !  two_bod_alpha_beta_mo_contracted(n,m,j,i) = \sum_{k,l} <ij|kl> <\Psi|a^{dagger}_{k,alpha} a^{dagger}_{l,beta} a_{n,beta} a_{m,alpha}|\Psi>
 END_DOC
 integer :: i,j,k,l,istate,m,n
 double precision :: wall_0,wall_1,accu
 two_bod_alpha_beta_mo_contracted = 0.d0
 print*,'providing two_bod_alpha_beta_mo_contracted...'
 call wall_time(wall_0)
 i=1 
 j=1
 double precision, allocatable :: integrals_array(:,:)
 allocate(integrals_array(mo_tot_num,mo_tot_num))
 call get_mo_bielec_integrals_ij(i,j,mo_tot_num,integrals_array,mo_integrals_map) 
 deallocate(integrals_array)

 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (l,k,i,j,n,m,accu,integrals_array,istate,mo_integrals_map) &
 !$OMP SHARED  (mo_tot_num,two_bod_alpha_beta_mo_transposed,N_states,two_bod_alpha_beta_mo_contracted) 
 allocate(integrals_array(mo_tot_num,mo_tot_num))
 do istate = 1, N_states 
 !$OMP DO              
  do i = 1, mo_tot_num ! 1 
   do j = 1, mo_tot_num  ! 2 
                                 !  1 2 
    call get_mo_bielec_integrals_ij(i,j,mo_tot_num,integrals_array,mo_integrals_map) 
    do m = 1, mo_tot_num ! 1 
     do n = 1, mo_tot_num ! 2 
      accu = 0.d0
      do l = 1, mo_tot_num ! 2 
       do k = 1, mo_tot_num ! 1 
        !                                        
        !                                        1 2 2 1                           1 2
        accu += two_bod_alpha_beta_mo_transposed(k,l,n,m,istate) * integrals_array(k,l) 
       enddo
      enddo
      !                                2 1 2 1           
      two_bod_alpha_beta_mo_contracted(n,m,j,i,istate) = accu 
     enddo
    enddo
   enddo
  enddo
 !$OMP END DO
 enddo
 deallocate(integrals_array)
 !$OMP END PARALLEL

 call wall_time(wall_1)
 print*,'two_bod_alpha_beta_mo_contracted provided in',dabs(wall_1-wall_0)

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, two_bod_alpha_beta_mo_transposed, (mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num,N_states)]
 implicit none
 BEGIN_DOC
 !  two_bod_alpha_beta_mo_transposed(i,j,k,l) = <Psi| a^{dagger}_{l,alpha} a^{dagger}_{k,beta} a_{j,beta} a_{i,alpha} | Psi>
 !                                   1 2 2 1
 !  note that no 1/2 factor is introduced in order to take into acccount for the spin symmetry
 END_DOC
 integer :: i,j,k,l,istate
 double precision :: cpu_0,cpu_1
 double precision :: accu
 two_bod_alpha_beta_mo_transposed= 0.d0
 print*,'providing two_bod_alpha_beta transposed ...'
 call cpu_time(cpu_0)
 do istate = 1, N_states 
  do i = 1, mo_tot_num
   do j = 1, mo_tot_num
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      !                                1 2 2 1                                 1 1 2 2
      two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) = two_bod_alpha_beta_mo(i,l,j,k,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
 call cpu_time(cpu_1)
 print*,'two_bod_alpha_beta transposed provided in',dabs(cpu_1-cpu_0)

 END_PROVIDER 

