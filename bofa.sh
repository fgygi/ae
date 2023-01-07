#!/bin/bash
# extract values of b vs a
for a in {1..50}
do src/vaepot 1 $a  2>/dev/null | grep b= |\
 awk '{print $2,$3}'|sed "s/a=//" | sed "s/b=//"
done

