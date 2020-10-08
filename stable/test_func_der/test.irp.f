!do i = 1, Npts
! delta_n_i ! variation de la densité
! d_e_d_n   ! derivée fonctionnelle
! int += delta_n_i * d_e_d_n * weight
!enddo

!function d_e_d_n(r)
!do k = 1, n_mo
! do l = 1, n_mo
!  d_e_d_n += v(l,k) * mo_k(r) * mo_l(r)
! enddo
!enddo
!end


