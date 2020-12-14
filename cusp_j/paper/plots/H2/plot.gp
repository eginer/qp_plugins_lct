
set terminal pdf enhanced
set key spacing 1.3                                                                                                       

set output "H2_fci.pdf"
set ylabel "Total energy for H_{2}"
set xlabel "interatomic distance (a.u.)"
set key right bottom
set xrange [1:3]
set yrange [-1.18:-1.12]
set key font "10,15"                                                                                                      


plot  "data_H2_e_fci" u 1:2 w lp lw 2 title "Exact", "" u 1:3 w lp lw 2 title "FCI AVDZ", "" u 1:4 w lp lw 2 title "FCI AVTZ", "" u 1:5 w lp lw 2 title "FCI AVQZ"


set terminal pdf enhanced
set output "H2_lda.pdf"
set ylabel "Total energy for H_{2}"
set xlabel "interatomic distance (a.u.)"
set key right bottom
set xrange [1:3]
set yrange [-1.18:-1.12]

plot  "data_H2_e_lda" u 1:2 w lp lw 2 title "Exact", "" u 1:3 w lp lw 2 title "{/Symbol m}_{UEG} AVDZ", "" u 1:4 w lp lw 2 title "{/Symbol m}_{UEG} AVTZ", "" u 1:5 w lp lw 2 title "{/Symbol m}_{UEG} AVQZ"


set terminal pdf enhanced
set output "H2_de_fci.pdf"
set ylabel "Energy (a.u.)" font ",15"
set xlabel "interatomic distance (a.u.)" font ",15"
set key right bottom
set xrange [1:3]
unset yrange 

plot  "data_H2_de_fci" u 1:3 w lp lw 2 title "FCI AVDZ", "" u 1:4 w lp lw 2 title "FCI AVTZ", "" u 1:5 w lp lw 2 title "FCI AVQZ"


set terminal pdf enhanced
set output "H2_de_lda.pdf"
set ylabel "Energy (a.u.)" font ",15"
set xlabel "interatomic distance (a.u.)" font ",15"
set key right bottom
set xrange [1:3]
unset yrange 

plot  "data_H2_de_lda" u 1:3 w lp lw 2 title "{/Symbol m}_{UEG} AVDZ", "" u 1:4 w lp lw 2 title "{/Symbol m}_{UEG} AVTZ", "" u 1:5 w lp lw 2 title "{/Symbol m}_{UEG} AVQZ"

set terminal pdf enhanced
set output "H2_de_fci_lda.pdf"
set ylabel "E_{exact} - E_X (m a.u.)" font ",15"
set xlabel "interatomic distance (a.u.)" font ",15"
set key right top
set xrange [1:3]
unset yrange 
set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
set label "R_{eq}" at 1.405,15.0


plot  "data_H2_de_fci" u 1:3 w lp lw 2 pt 1 title "X: FCI AVDZ", "data_H2_de_lda" u 1:3 w lp lw 2 pt 4 title "X: {/Symbol m}_{UEG} AVDZ",\
      "data_H2_de_fci" u 1:4 w lp lw 2 pt 1 title "X: FCI AVTZ", "data_H2_de_lda" u 1:4 w lp lw 2 pt 4 title "X: {/Symbol m}_{UEG} AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"

set terminal pdf enhanced
set output "H2_de_fci_rsc.pdf"
set ylabel "E_{exact} - E_X (m a.u.)" font ",15"
set xlabel "interatomic distance (a.u.)" font ",15"
set key right top
set xrange [1:3]
unset yrange 
#set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
#set label "R_{eq}" at 1.405,15.0


plot  "data_H2_de_fci" u 1:3 w lp lw 2 pt 1 title "X: FCI AVDZ", "data_H2_de_rsc" u 1:3 w lp lw 2 pt 4 title "X: {/Symbol m}_{r_{s,c}} AVDZ",\
      "data_H2_de_fci" u 1:4 w lp lw 2 pt 1 title "X: FCI AVTZ", "data_H2_de_rsc" u 1:4 w lp lw 2 pt 4 title "X: {/Symbol m}_{r_{s,c}} AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"

set terminal pdf enhanced
set output "H2_de_fci_frogg.pdf"
set ylabel "E_{exact} - E_X (m a.u.)" font ",15"
set xlabel "interatomic distance (a.u.)" font ",15"
set key right top
set xrange [1:3]
unset yrange 
#set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
#set label "R_{eq}" at 1.405,15.0


plot  "data_H2_de_fci" u 1:3 w lp lw 2 pt 1 title "X: FCI AVDZ", "data_H2_de_ten_no" u 1:3 w lp lw 2 pt 4 title "X: FROGG AVDZ",\
      "data_H2_de_fci" u 1:4 w lp lw 2 pt 1 title "X: FCI AVTZ", "data_H2_de_ten_no" u 1:4 w lp lw 2 pt 4 title "X: FROGG AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"

set terminal pdf enhanced
set output "H2_de_fci_mu_ten_no.pdf"
set ylabel "E_{exact} - E_X (m a.u.)" font ",15"
set xlabel "interatomic distance (a.u.)" font ",15"
set key right top
set xrange [1:3]
unset yrange 
#set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
#set label "R_{eq}" at 1.405,15.0


plot  "data_H2_de_fci" u 1:3 w lp lw 2 pt 1 title "X: FCI AVDZ", "data_H2_de_mu_ten_no" u 1:3 w lp lw 2 pt 4 title "X: {/Symbol m}_{FROGG} AVDZ",\
      "data_H2_de_fci" u 1:4 w lp lw 2 pt 1 title "X: FCI AVTZ", "data_H2_de_mu_ten_no" u 1:4 w lp lw 2 pt 4 title "X: {/Symbol m}_{FROGG} AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"
