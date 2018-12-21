====================
Integrals_ijkl_in_r3
====================


This modules contains the following integrals: 
 \int d^3r i(r) j(r) k(r) l(r) 

on the |AO| and |MO| basis. This can be useful for some multi-determinant correlation energy in |RSDFT|. 

To fetch an |AO| integral, use the
`get_ao_bielec_integral_ijkl_r3(i,j,k,l,ao_integrals_ijkl_r3_map)` function, and
to fetch an |MO| integral, use
`get_mo_bielec_integral_ijkl_r3(i,j,k,l,mo_integrals_ijkl_r3_map)` 

As they are all symmetric, it does not matter the conventions. 
