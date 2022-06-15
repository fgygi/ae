#!/bin/bash
for f in $*
do
  awk '/<radial_function>/,/<\/radial_function>/' $f | \
  awk 'NR>2 {print l;} {l=$0}' >> vaeplot$$
  echo " " >> vaeplot$$
  echo " " >> vaeplot$$
done
gnuplot -persist <<EOF
set grid
plot "vaeplot$$" w l
EOF
rm vaeplot$$
