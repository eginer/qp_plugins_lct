source /home/emmanuel/qp2/quantum_package.rc
qp create_ezfio -b aug-cc-pvdz He.xyz 
qp run scf 
qp set determinants threshold_save_wf 0.
qp run cis 
qp run diag_dress_iter | tee He.cis_tc_0.out
for i in `seq 1 10`
 do
 qp run tc_scf | tee He.ezfio.tc_scf_${i}.out
 qp run diag_dress_iter | tee He.cis_tc_${i}.out
done
