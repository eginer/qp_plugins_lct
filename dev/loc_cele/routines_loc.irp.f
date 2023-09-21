subroutine loc_cele_routine_general(n_orb_for_loc,index_orb_for_loc,mo_guess_for_loc)
     implicit none
     integer, intent(in) :: n_orb_for_loc,index_orb_for_loc(n_orb_for_loc)
     double precision, intent(in) :: mo_guess_for_loc(ao_num,n_orb_for_loc)
     integer id1,i_atom,shift,shift_h
     character*1 jobz,uplo
     character*64 file1,file2
     character*72 cdum
     integer n,i,j,k,l,isym,nsym,idum,m
     integer info,lwork
     logical *1 z54
     integer                        :: nmo0(8), nrot(8), index_rot(1000,1)
     double precision               :: accu_norm

     double precision, allocatable  :: cmo(:,:,:),cmoref(:,:,:),newcmo(:,:,:)
     double precision, allocatable  :: s(:,:,:),dum,ddum(:,:),ovl(:,:)
     double precision, allocatable  :: w(:),work(:),t(:,:),wi(:,:)
     integer, allocatable           :: irot(:,:),ipiv(:)
     character*72, allocatable      :: string(:,:)
     id1 = 700

     allocate( cmo(id1,id1,1),cmoref(id1,id1,1),newcmo(id1,id1,1))
     allocate( s(id1,id1,1),dum,ddum(id1,id1),ovl(id1,id1))
     allocate( w(id1),work(3*id1),t(id1,id1),wi(id1,id1), irot(id1,8),string(id1,8),ipiv(id1))
     print*,'coucou'
     z54=.false.
     accu_norm = 0.d0
     do i =1,mo_num
       accu_norm += dabs(mo_overlap(i,i))
     enddo
     print*,'accu_norm = ',accu_norm
     nsym = 1
     nmo0(1) = mo_num
     cmo = 0.d0
     do isym=1,nsym
       do i=1,nmo0(isym)
         do j = 1, ao_num
           cmo(j,i,isym) = mo_coef(j,i)
         enddo
       enddo
     enddo
     do isym=1,nsym
       do j=1,mo_num
         do i=1,ao_num
           newcmo(i,j,isym)=cmo(i,j,isym)
         enddo
       enddo
     enddo
     nrot(1) = n_orb_for_loc ! number of orbitals to be localized
     
     
     
     
     cmoref = 0.d0
     irot = 0
     
  do i = 1, n_orb_for_loc
   irot(i,1) = index_orb_for_loc(i)
  enddo
  do i = 1, n_orb_for_loc
   do j = 1, ao_num
    cmoref(j,i,1)   = mo_guess_for_loc(j,i)
   enddo
  enddo

  do i = 1, nrot(1)
    print*,'irot(i,1) = ',irot(i,1)
  enddo
  
     
     
     do i = 1, ao_num
       do j =1, ao_num
         s(j,i,1) =  ao_overlap(j,i)
       enddo
     enddo
     !Now big loop over symmetry
     
     
     
     do isym=1,nsym
       print*,'isym = ',isym 
       if (nrot(isym).eq.0) cycle
       do i=1,nrot(isym) ! mos
         do j=1,nrot(isym)! ref vector 
           ovl(j,i)=0.d0
           print*,i,j
            do m = 1, ao_num
             do k = 1, ao_num
               ovl(j,i)=ovl(j,i)+cmoref(k,j,isym)*cmo(m,irot(i,isym),isym) * s(k,m,isym)
             enddo
            enddo
         enddo
       enddo
       print*,'OVL === '
       do i = 1, nrot(isym)
        write(*,*)ovl(1:nrot(isym),i)
       enddo
       call maxovl(nrot(isym),nrot(isym),ovl,t,wi)
       do i=1,nrot(isym)
         do j=1,ao_num
           !         write (6,*) 'isym,',isym,nrot(isym),nmo0(isym)
           newcmo(j,irot(i,isym),isym)=0.d0
           do k=1,nrot(isym)
             newcmo(j,irot(i,isym),isym)=newcmo(j,irot(i,isym),isym) + cmo(j,irot(k,isym),isym)*t(k,i)
           enddo
         enddo
       enddo
       
     enddo !big loop over symmetry
     
     10 format (4E19.12)
     
     
     !  Now we copy the newcmo into the mo_coef
     
     mo_coef = 0.d0
     do isym=1,nsym
       do i=1,nmo0(isym)
         do j = 1, ao_num
           mo_coef(j,i) = newcmo(j,i,isym)
         enddo
       enddo
     enddo
     !      pause
     
     
     ! we say that it hase been touched, and valid and that everything that
     ! depends on mo_coef must not be reprovided
     accu_norm = 0.d0
     do i =1,mo_num
       accu_norm += dabs(mo_overlap(i,i))
     enddo
     print*, 'accu_norm = ',accu_norm
     ! We call the routine that saves mo_coef in the ezfio format
     call save_mos
     
end

