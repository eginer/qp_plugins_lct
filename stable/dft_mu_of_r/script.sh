ezfio=$1
ezfio=${ezfio%/}
nstate=`cat $ezfio/determinants/n_states`

for i in `seq 1 $nstate`
 do
 echo $i 
 state=${ezfio}_state_$i
 cp -r $ezfio $state 
 qp_edit -s $i $state 
#qp_run truncate_wf_spin_no_H $state 
 qp_run all_mu_of_r_corrections_md $state | tee ${state}.DFT.out 
done
