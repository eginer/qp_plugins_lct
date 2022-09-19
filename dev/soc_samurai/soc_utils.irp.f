integer function how_many_sp(S,ms)
 implicit none
 BEGIN_DOC
 ! routine that returns how many S^+ one can act on a state |S,ms>
 END_DOC
 double precision, intent(in) :: S,ms
 double precision :: ms_tmp,factor_s_p,coef
 integer :: itmp 
 itmp = 0
 coef = 1.d0
 ms_tmp = ms
 do while (dabs(coef).gt.1.d-10)
  coef = factor_s_p(S,ms_tmp)
  if(dabs(coef).gt.1.d-10)then
   itmp += 1
  endif
  ms_tmp += 1.d0
 enddo
 how_many_sp = itmp
end


integer function how_many_sm(S,ms)
 implicit none
 BEGIN_DOC
 ! routine that returns how many S^- one can act on a state |S,ms>
 END_DOC
 double precision, intent(in) :: S,ms
 double precision :: ms_tmp,factor_s_m,coef
 integer :: itmp 
 itmp = 0
 coef = 1.d0
 ms_tmp = ms
 do while (dabs(coef).gt.1.d-10)
  coef = factor_s_m(S,ms_tmp)
  if(dabs(coef).gt.1.d-10)then
   itmp += 1
  endif
  ms_tmp -= 1.d0
 enddo
 how_many_sm = itmp
end

double precision function factor_s_p(S,ms)
 implicit none
 double precision, intent(in) :: S,ms
 factor_s_p = dsqrt(S*(S+1.d0) - ms*(ms+1.d0))
end

double precision function factor_s_m(S,ms)
 implicit none
 double precision, intent(in) :: S,ms
 factor_s_m = dsqrt(S*(S+1.d0) - ms*(ms-1.d0))
end
