#!/bin/bash
# extract values of b vs a
for a in {36..180}
do src/vaepot_1.0 18 $a  2>/dev/null | grep b= |\
 awk '{print $2,$3}'|sed "s/a=//" | sed "s/b=//"
done

