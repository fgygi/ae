#!/bin/bash
# solve the AE problem for Z in a range [amin,amax]
for a in $(seq 1 1 8)
  do
  src/aecheck 1 $a 0 100 0.0001
done

