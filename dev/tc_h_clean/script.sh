 source /home/emmanuel/qp2/quantum_package.rc
 atom=He
 qp run scf 
 echo $atom > ${atom}.xyz
 qp create_ezfio -b aug-cc-pvdz ${atom}.xyz 
 qp run scf 
 # prints out typical values of mu, mu_lda and mu_rsc_lda seem good 
 qp run print_mu_av_tc | tee ${EZFIO_FILE}.mu_av.out 
 mu_rsc_lda=` grep "average_mu_rs_c_lda" ${EZFIO_FILE}.mu_av.out | cut -d "=" -f 2`
 qp set ao_two_e_erf_ints mu_erf $mu_rsc_lda 
 # you set the threshold to ZERO to save all the determinants in the CIS 
 qp set determinants threshold_save_wf 0.
 # You diagonalize the TC Hamiltonian within the set of determinants saved in the EZFIO
for i in `seq 1 10`
 do
 qp run cis # list of determinants CIS 
 qp run diag_dress_iter | tee ${EZFIO_FILE}.tc_cis_${i}.out # diagonalize Htc in the basis CIS -> EZFIO 
 # You save the natural orbitals right eigenvector of CIS-TC
 qp run save_natorb | tee ${EZFIO_FILE}.tc_natorb_${i}.out 
 # you iterate ...
done
