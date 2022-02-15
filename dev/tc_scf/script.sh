 source /home/emmanuel/qp2/quantum_package.rc
 atom=Be
 echo $atom > ${atom}.xyz
 rm -rf $atom.ezfio
# qp create_ezfio -b "O:cc-pcvdz|H:cc-pvdz" ${atom}.xyz 
 qp create_ezfio -b cc-pcvdz ${atom}.xyz -m 1 -o 
# qp create_ezfio -b basis ${atom}.xyz 
 qp run scf 
 # prints out typical values of mu, mu_lda and mu_rsc_lda seem good 
 qp run print_mu_av_tc | tee ${EZFIO_FILE}/mu_av.out 
 mu_rsc_lda=` grep "average_mu_rs_c_lda" ${EZFIO_FILE}/mu_av.out | cut -d "=" -f 2`
# qp set tc_h_clean full_tc_h_solver True
 qp set ao_two_e_erf_ints mu_erf $mu_rsc_lda 
 # you set the threshold to ZERO to save all the determinants in the CIS 
 qp set determinants save_threshold 0.
 qp set tc_h_clean three_body_h_tc true 
 qp set tc_h_clean pure_three_body_h_tc False
 qp set tc_h_clean double_normal_ord False
 cp -r $EZFIO_FILE ${EZFIO_TILE}_tc_scf
 qp set_file ${EZFIO_TILE}_tc_scf 
 qp run tc_scf | tee ${EZFIO_FILE}.out 
 qp set_file $atom.ezfio 

# qp run cisd 
# qp run diag_dress_iter | tee ${EZFIO_FILE}/tc_cisd_0.out # diagonalize Htc in the basis CIS -> EZFIO 

 # You diagonalize the TC Hamiltonian within the set of determinants saved in the EZFIO
for i in `seq 1 10`
 do
 qp run cis # list of determinants CIS 
 qp run diag_dress_iter | tee ${EZFIO_FILE}/tc_cis_${i}.out # diagonalize Htc in the basis CIS -> EZFIO 
 qp run print_wf | tee ${EZFIO_FILE}/tc_cis_wf_${i}.out 
 # You save the natural orbitals right eigenvector of CIS-TC
 qp run save_natorb | tee ${EZFIO_FILE}/tc_natorb_${i}.out 
 # you iterate ...
# qp run cisd 
# qp run diag_dress_iter | tee ${EZFIO_FILE}/tc_cisd_${i}.out # diagonalize Htc in the basis CIS -> EZFIO 
# qp run print_wf | tee ${EZFIO_FILE}/tc_cisd_wf_${i}.out 
done
# qp run cisd 
# qp run diag_dress_iter | tee ${EZFIO_FILE}/tc_cisd_${i}.out # diagonalize Htc in the basis CIS -> EZFIO 
