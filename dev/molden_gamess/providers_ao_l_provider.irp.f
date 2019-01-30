 BEGIN_PROVIDER [ integer, ao_l_powers, (ao_num,10) ]
 implicit none
 integer :: i
 character*(4) :: give_ao_character_space
 do i=1,ao_num

  if(ao_l(i)==0)then
  ! S type AO
  elseif(ao_l(i) == 1)then
  ! P type AO
   if(ao_power(i,1)==1)then
    ao_l_powers(i,1) = 1
   elseif(ao_power(i,2) == 1)then
    ao_l_powers(i,1) = 2
   else
    ao_l_powers(i,1) = 3
   endif
  elseif(ao_l(i) == 2)then
  ! D type AO
   if(ao_power(i,1)==2)then
    ao_l_powers(i,1) = 1
    ao_l_powers(i,2) = 1
   elseif(ao_power(i,2) == 2)then
    ao_l_powers(i,1) = 2
    ao_l_powers(i,2) = 2
   elseif(ao_power(i,3) == 2)then
    ao_l_powers(i,1) = 3
    ao_l_powers(i,2) = 3
   elseif(ao_power(i,1) == 1 .and. ao_power(i,2) == 1)then
    ao_l_powers(i,1) = 1
    ao_l_powers(i,2) = 2
   elseif(ao_power(i,1) == 1 .and. ao_power(i,3) == 1)then
    ao_l_powers(i,1) = 1
    ao_l_powers(i,2) = 3
   else
    ao_l_powers(i,1) = 2
    ao_l_powers(i,2) = 3
   endif
  elseif(ao_l(i) == 3)then
  ! F type AO
   if(ao_power(i,1)==3)then
    ao_l_powers(i,1) = 1
    ao_l_powers(i,2) = 1
    ao_l_powers(i,3) = 1
   elseif(ao_power(i,2) == 3)then
    ao_l_powers(i,1) = 2
    ao_l_powers(i,2) = 2
    ao_l_powers(i,3) = 2
   elseif(ao_power(i,3) == 3)then
    ao_l_powers(i,1) = 3
    ao_l_powers(i,2) = 3
    ao_l_powers(i,3) = 3
   elseif(ao_power(i,1) == 2 .and. ao_power(i,2) == 1)then
    ao_l_powers(i,1) = 1
    ao_l_powers(i,2) = 1
    ao_l_powers(i,3) = 2
   elseif(ao_power(i,1) == 2 .and. ao_power(i,3) == 1)then
    ao_l_powers(i,1) = 1
    ao_l_powers(i,2) = 1
    ao_l_powers(i,3) = 3
   elseif(ao_power(i,2) == 2 .and. ao_power(i,1) == 1)then
    ao_l_powers(i,1) = 2
    ao_l_powers(i,2) = 2
    ao_l_powers(i,3) = 1
   elseif(ao_power(i,2) == 2 .and. ao_power(i,3) == 1)then
    ao_l_powers(i,1) = 2
    ao_l_powers(i,2) = 2
    ao_l_powers(i,3) = 3
   elseif(ao_power(i,3) == 2 .and. ao_power(i,1) == 1)then
    ao_l_powers(i,1) = 3
    ao_l_powers(i,2) = 3
    ao_l_powers(i,3) = 1
   elseif(ao_power(i,3) == 2 .and. ao_power(i,2) == 1)then
    ao_l_powers(i,1) = 3
    ao_l_powers(i,2) = 3
    ao_l_powers(i,3) = 2
   elseif(ao_power(i,3) == 1 .and. ao_power(i,2) == 1 .and. ao_power(i,3) == 1)then
    ao_l_powers(i,1) = 1
    ao_l_powers(i,2) = 2
    ao_l_powers(i,3) = 3
   endif
  elseif(ao_l(i) .ge. 4)then
   print*,'Warning !! the correct ordering of the AOs for the GAMESS(US) format is not done' 
   print*,'for angular momentum higer than F !!'
   stop
  endif
 enddo
END_PROVIDER

