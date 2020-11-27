
set key spacing 1.3                                                                                                       
set terminal pdf enhanced
set output "H2_de_fci_mu_0_3.pdf"
set ylabel "E_{exact} - E_X (m a.u.)"
set xlabel "interatomic distance (a.u.)"
set key right top
set xrange [1:3]
unset yrange 
set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
set label "R_{eq}" at 1.401,15.0


plot  "data_H2_de_fci" u 1:3 w lp pt 1 title "X: FCI AVDZ", "data_H2_de_mu_0.3" u 1:3 w lp pt 4 title "X: {/Symbol m} = 0.3 AVDZ",\
      "data_H2_de_fci" u 1:4 w lp pt 1 title "X: FCI AVTZ", "data_H2_de_mu_0.3" u 1:4 w lp pt 4 title "X: {/Symbol m} = 0.3 AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"


set key spacing 1.3                                                                                                       
set terminal pdf enhanced
set output "H2_de_fci_mu_0_5.pdf"
set ylabel "E_{exact} - E_X (m a.u.)"
set xlabel "interatomic distance (a.u.)"
set key right top
set xrange [1:3]
unset yrange 
set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
set label "R_{eq}" at 1.401,15.0


plot  "data_H2_de_fci" u 1:3 w lp pt 1 title "X: FCI AVDZ", "data_H2_de_mu_0.5" u 1:3 w lp pt 4 title "X: {/Symbol m} = 0.5 AVDZ",\
      "data_H2_de_fci" u 1:4 w lp pt 1 title "X: FCI AVTZ", "data_H2_de_mu_0.5" u 1:4 w lp pt 4 title "X: {/Symbol m} = 0.5 AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"

set key spacing 1.3                                                                                                       
set terminal pdf enhanced
set output "H2_de_fci_mu_0_7.pdf"
set ylabel "E_{exact} - E_X (m a.u.)"
set xlabel "interatomic distance (a.u.)"
set key right top
set xrange [1:3]
unset yrange 
set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
set label "R_{eq}" at 1.401,15.0


plot  "data_H2_de_fci" u 1:3 w lp pt 1 title "X: FCI AVDZ", "data_H2_de_mu_0.7" u 1:3 w lp pt 4 title "X: {/Symbol m} = 0.7 AVDZ",\
      "data_H2_de_fci" u 1:4 w lp pt 1 title "X: FCI AVTZ", "data_H2_de_mu_0.7" u 1:4 w lp pt 4 title "X: {/Symbol m} = 0.7 AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"

set key spacing 1.3                                                                                                       
set terminal pdf enhanced
set output "H2_de_fci_mu_0_9.pdf"
set ylabel "E_{exact} - E_X (m a.u.)"
set xlabel "interatomic distance (a.u.)"
set key right top
set xrange [1:3]
unset yrange 
set arrow from 1.401, graph 0 to 1.401, graph 1 nohead
set label "R_{eq}" at 1.401,15.0


plot  "data_H2_de_fci" u 1:3 w lp pt 1 title "X: FCI AVDZ", "data_H2_de_mu_0.9" u 1:3 w lp pt 4 title "X: {/Symbol m} = 0.9 AVDZ",\
      "data_H2_de_fci" u 1:4 w lp pt 1 title "X: FCI AVTZ", "data_H2_de_mu_0.9" u 1:4 w lp pt 4 title "X: {/Symbol m} = 0.9 AVTZ",\
      -1.6 dt 4 lt 6 notitle, 1.6 dt 4 lt 6 title "chemical accuracy"
