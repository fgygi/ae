#!/bin/bash
# get error from aecheck output
# use: ./get_err.sh file.dat
for n in E_2 E_3 E_4 E_5
do
  echo '# ' $n
  paste <( grep 'a = ' $1 | awk '{print $6}' )  \
        <( grep $n $1 | awk '{print $5}' )
  echo; echo
done

