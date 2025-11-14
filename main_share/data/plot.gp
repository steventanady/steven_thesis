if (ARG1 eq "pdf") {
  set terminal push
  set terminal pdfcairo enhanced colour dashed \
    font "Liberation Serif,12" linewidth 2
  set out "data.pdf"
}
offset_1 = system ("grep offset pantheonSample.d | sed 's/.*= *//' | awk '{ printf (\"%5.2f\", $0) }'")
offset_2 = system ("grep offset baoSample.d | sed 's/.*= *//' | awk '{ printf (\"%5.2f\", $0) }'")
set xlabel "redshift z"
set ylabel "distance modulus"
set key bottom right
set grid
plot \
  "unionSample2.d" u 1:2:($2+$3):($2-$4) \
  ls 1 w yerr t "Union 2", \
  "pantheonSample.d" u 1:($2+offset_1):($2+offset_1+$3):($2+offset_1-$4) \
  ls 2 w yerr t "Pantheon", \
  "baoSample.d" u 1:($2+offset_2):($2+offset_2+$3):($2+offset_2-$4) \
  ls 3 w yerr t "BAO"
if (ARG1 eq "pdf") {
  set out
  set terminal pop
}
