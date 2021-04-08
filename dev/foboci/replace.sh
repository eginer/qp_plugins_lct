
word='mo_tot_num'
neword='mo_num'

word='mo_mono_elec_integral'
neword='mo_one_e_integrals'

word='one_body_dm_mo_alpha_average'
neword='one_e_dm_mo_alpha_average'

word='one_body_dm_mo_beta_average'
neword='one_e_dm_mo_beta_average'

word='one_body'
neword='one_e'

word='do_mono_excitation'
neword='do_single_excitation'

word='davidson_criterion'
neword='threshold_davidson'

word='mono'
neword='single'

list=`grep -i $word *.f | cut -d ":" -f 1 | uniq`
for i in $list
 do 
  echo $i
  sed -e "s/`echo $word`/`echo $neword`/g" $i > pouet 
  mv pouet $i 
done
