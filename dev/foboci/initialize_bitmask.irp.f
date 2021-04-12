use bitmasks

subroutine initialize_bitmask_to_restart_ones                                                                                                
 implicit none
 integer :: i,j,l,m
 integer :: ispin
  BEGIN_DOC
 ! Initialization of the generators_bitmask to the restart bitmask
  END_DOC
 do i = 1, N_int
   do ispin=1,2
     generators_bitmask(i,ispin,s_hole ) =   generators_bitmask_restart(i,ispin,s_hole )
     generators_bitmask(i,ispin,s_part ) =   generators_bitmask_restart(i,ispin,s_part )
     generators_bitmask(i,ispin,d_hole1) =   generators_bitmask_restart(i,ispin,d_hole1)
     generators_bitmask(i,ispin,d_part1) =   generators_bitmask_restart(i,ispin,d_part1)
     generators_bitmask(i,ispin,d_hole2) =   generators_bitmask_restart(i,ispin,d_hole2)
     generators_bitmask(i,ispin,d_part2) =   generators_bitmask_restart(i,ispin,d_part2)
   enddo
 enddo
end



BEGIN_PROVIDER [ integer(bit_kind), generators_bitmask_restart, (N_int,2,6) ]
 implicit none
 BEGIN_DOC
 ! Bitmasks for generator determinants.
 ! (N_int, alpha/beta, hole/particle, generator).
 !
 ! 3rd index is :
 !
 ! * 1 : hole     for single exc
 !
 ! * 2 : particle for single exc
 !
 ! * 3 : hole     for 1st exc of double
 !
 ! * 4 : particle for 1st exc of double
 !
 ! * 5 : hole     for 2nd exc of double
 !
 ! * 6 : particle for 2nd exc of double
 !
 END_DOC
 logical                        :: exists
 PROVIDE N_int

 integer :: i,ispin
 integer :: j
 integer, save :: ifirst = 0
 if(ifirst == 0)then
  ifirst = 1
  do ispin=1,2
      do i=1,N_int
        generators_bitmask_restart(i,ispin,s_hole ) = reunion_of_inact_act_bitmask(i,ispin)
        generators_bitmask_restart(i,ispin,s_part ) = reunion_of_act_virt_bitmask(i,ispin)
        generators_bitmask_restart(i,ispin,d_hole1) = reunion_of_inact_act_bitmask(i,ispin)
        generators_bitmask_restart(i,ispin,d_part1) = reunion_of_act_virt_bitmask(i,ispin)
        generators_bitmask_restart(i,ispin,d_hole2) = reunion_of_inact_act_bitmask(i,ispin)
        generators_bitmask_restart(i,ispin,d_part2) = reunion_of_act_virt_bitmask(i,ispin)
      enddo
  enddo
 else
  print*,'PB in generators_bitmask_restart !!!'
 endif

END_PROVIDER

