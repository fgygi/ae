#!/bin/bash
# solve the AE problem for Z in a range [amin,amax]
for a in $(seq 24 4 48)
  do
  src/aecheck 6 $a 0 30 0.0001
done

