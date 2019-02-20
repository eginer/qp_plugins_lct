#!/usr/bin/env qpsh
 ezfio=$1
 qp set_file $ezfio

 x=$2
 y=$3
 z=$4
 delta_r12=$5

 qp set firth_order_der x_center_wee $x
 qp set firth_order_der y_center_wee $y
 qp set firth_order_der z_center_wee $z
 qp set firth_order_der delta_r12    $delta_r12

 qp run print_wee > pouet
 out=${ezfio}_${x}_${y}_${z}_dr12_$delta_r12
 cp fort.33 $out
 mu=`grep "mu_ana =" pouet | cut -d "=" -f 2`


cat << EOF > ${out}.plt
set term pdf
set output "$out.pdf"
set xlabel "r_{12}"
set ylabel "Interaction"
ymax=2*$mu/sqrt(pi) + 1.
set yrange [0:ymax]

plot "${out}" w l title "W_{ee} paper", "" u 1:3 w l title "Fit mu_{num}" , "" u 1:4 w l title "Fit mu_{ana}", 1/x title "exact"

EOF

gnuplot ${out}.plt
