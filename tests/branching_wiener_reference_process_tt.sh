#!/bin/sh

for fn
do
  grep "#M_areatt " $fn | awk ' BEGIN {max_slot=-1; } { slot=$3; if (slot>max_slot) { max_slot=slot; } tt[slot]=$4; ell=$5; hopping=$6; epsilon=$7; sigma=$8; dim=$9; for (i=0; i<=8; i++) {v=$(i+14); m1[slot,i]+=v; m2[slot,i]+=(v*v);} } END { for (slot=0; slot<=max_slot; slot++) { printf("%i %g %i %i %g %g %g %i", slot, tt[slot], m1[slot,0], ell, hopping, epsilon, sigma, dim);  for (i=1; i<=8; i++) { m1[slot,i]/=m1[slot,0]; m2[slot,i]/=m1[slot,0]; printf(" %g %g", m1[slot,i], sqrt((m2[slot,i]-m1[slot,i]*m1[slot,i])/(m1[slot,0]-1))); } printf("\n");}} '
done


