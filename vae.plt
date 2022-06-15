#!/bin/bash
for f in $*
do
  awk '/<radial_potential>/,/<\/radial_potential>/' $f | \
  awk 'NR>2 {print l;} {l=$0}' >> vaeplot$$
  echo " " >> vaeplot$$
  echo " " >> vaeplot$$
done
gnuplot -persist <<EOF
set grid
plot "vaeplot$$" w l
EOF
rm vaeplot$$
