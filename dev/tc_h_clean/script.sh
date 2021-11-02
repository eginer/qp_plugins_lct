 qp run scf 
 qp create_ezfio -b cc-pcvdz Ne.xyz 
 qp run scf 
 # prints out typical values of mu, mu_lda and mu_rsc_lda seem good 
 qp run print_mu_av_tc 
 qp set ao_two_e_erf_ints mu_erf 1.5534429522451609 
 # you set the threshold to ZERO to save all the determinants in the CIS 
 qp set determinants threshold_save_wf 0.
 qp run cis 
 # You diagonalize the TC Hamiltonian within the set of determinants saved in the EZFIO
 qp run diag_dress_iter | tee Ne.ezfio.tc_cis.out 
 # You save the natural orbitals 
 qp run save_natorb 
 # you iterate ...
 qp run cis 
 qp run diag_dress_iter | tee Ne.ezfio.tc_cis_2.out 

