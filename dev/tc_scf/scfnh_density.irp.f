BEGIN_PROVIDER [double precision, Q_ao, (ao_num,ao_num) ]

  BEGIN_DOC
  ! eq (17) in P Salvador, Int J Quantum Chem, 2009
  END_DOC

  call dgemm( 'N', 'T', ao_num, ao_num, mo_num, 2.d0  &
             , Reig_tcFock_mo, size(Reig_tcFock_mo,1) &
             , Leig_tcFock_mo, size(Leig_tcFock_mo,1) &
             , 0.d0, Q_ao, size(Q_ao,1) )

END_PROVIDER

