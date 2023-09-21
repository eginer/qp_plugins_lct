! -*- F90 -*-
 BEGIN_PROVIDER [integer, Aindx_mpq, (n_det,n_act_orb,n_act_orb,2)]
&BEGIN_PROVIDER [real*8, Aval_mpq, (n_det,n_act_orb,n_act_orb,2)]
&BEGIN_PROVIDER [logical, secnd_Aval_mpq, (n_det,n_act_orb,n_act_orb)]
BEGIN_DOC
! formula tape for 1-electron integrals
! 
! we have to apply E_pq = destroy q, create p for spin alpha and spin
! beta  E_pq = (a+_p a_q) + (a+_-p a_-q)
! entries are 2, 1, 0, -1
! very few non-zero entries
!
! call do_spinfree_mono_excitation does the job

! a+_i a_j |n> = (1-n_i+delta_ij) n_j e_ij Gamma(n)_iGamma(n)_j|k>
!
! k_l=n_l if l neq i,j
! k_i = 1, k_j=delta_ij, e_ij = -1 for i>j; 1 oherwise
!
! we generate directly the packed form 
! Amnpq(mu,nu,p,q) -> Aindx_mpq(mu,p,q,1)=nu and the value goes to Aval_mpq
!    if there is a second one, then secnd_Aval_mpq is true
! as p.ne.q can generate different determinants for alpha and beta strings if
! p and q are singly occupied in determinant mu
!
END_DOC
      implicit none
      integer :: i,j,k,l,p,q,mu,nu,ii
      integer :: ierr,ierr1,ierr2,nu1,nu2,ihole,ipart
      integer(bit_kind), allocatable :: det_mu(:,:)
      integer(bit_kind), allocatable :: det_mu_ex1(:,:)
      integer(bit_kind), allocatable :: det_mu_ex2(:,:)
 
      real*8 :: phase,phase1,phase2
      allocate(det_mu(N_int,2))
      allocate(det_mu_ex1(N_int,2))
      allocate(det_mu_ex2(N_int,2))

      write(6,*) ' creating the packed formula tape Amnpq '
      write(6,*) ' n_det = ',n_det

      Aindx_mpq=-1
      Aval_mpq=0.D0
      secnd_Aval_mpq=.false.

integer :: n_one
      n_one=0
        do mu=1,n_det
         det_mu(1:N_int,1:2)=psi_det(1:N_int,1:2,mu)
!        call det_extract(det_mu,mu,N_int)
!         write(6,*) ' Determinant No ',mu
!         call print_det(det_mu,N_int)
        end do

!       do mu=1,n_det
!       write(6,*) 'F-tape det No ',mu,psi_coef(mu,1),Dkpq(mu,3,1)
!       det_mu(1:N_int,1:2)=psi_det(1:N_int,1:2,mu)
!       call det_extract(det_mu,mu,N_int)
!       call print_det(det_mu,N_int)
!       end do
! apply E_pq
       do p=1,n_act_orb
        ipart=list_act(p)
        do q=1,n_act_orb
         ihole=list_act(q)
         do mu=1,n_det
          det_mu(1:N_int,1:2)=psi_det(1:N_int,1:2,mu)
          call det_extract(det_mu,mu,N_int)
!         call det_copy(det_mu,det_mu_ex1,N_int)
!         call det_copy(det_mu,det_mu_ex2,N_int)
          call do_spinfree_mono_excitation(det_mu,det_mu_ex1  &
              ,det_mu_ex2,nu1,nu2,ihole,ipart,phase1,phase2,ierr1,ierr2)
          if (ierr1.eq.1) then
           if (nu1.ne.-1) then
            Aindx_mpq(mu,p,q,1)=nu1
            Aval_mpq(mu,p,q,1)+=phase1
            n_one+=1
!          else
!           write(6,*) ' I excitation not possible ? strange ... ',mu,p,q
!           write(6,*) ' ihole, ipart : ',ihole,ipart
!           call print_det(det_mu,N_int)
!           stop ' excitation not possible ? strange ... '
           end if
          end if
          if (ierr2.eq.1) then
           if (nu2.ne.-1) then
            if ((nu1.eq.nu2).or.(nu1.eq.-1)) then
             Aindx_mpq(mu,p,q,1)=nu2
             Aval_mpq(mu,p,q,1)+=phase2
             secnd_Aval_mpq(mu,p,q)=.false.
             if (nu1.eq.-1) then
              n_one+=1
             end if
            else
             Aindx_mpq(mu,p,q,2)=nu2
             Aval_mpq(mu,p,q,2)+=phase2
             secnd_Aval_mpq(mu,p,q)=.true.
             n_one+=1
            end if
!          else
!           write(6,*) ' II excitation not possible ? strange ... ',mu,p,q
!           stop ' excitation not possible ? strange ... '
           end if
          end if
         end do
        end do
       end do
      
      write(6,*) ' filled 1e formula tape with ',n_one,' non-zero entries'
!
!     open(unit=12,file='Apacked_mpq.dat',form='formatted',status='unknown')
!     do mu=1,n_det
!      do p=1,n_act_orb
!       do q=1,n_act_orb
!        do i=1,2
!         if (Aindx_mpq(mu,p,q,i).ne.-1) then
!             write(12,'(i2,I8,2I4,i8,F6.1)') i,mu,p,q   &
!                 ,Aindx_mpq(i,mu,p,q),Aval_mpq(mu,p,q,i)
!         end if
!        end do
!       end do
!      end do
!     end do
!     close(12)

END_PROVIDER

