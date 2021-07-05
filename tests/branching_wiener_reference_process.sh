#!/bin/sh

for fn
do
  grep "#M_area " $fn | awk ' { ell=$3; hopping=$4; epsilon=$5; sigma=$6; dim=$7; for (i=0; i<=8; i++) {v=$(i+12); m1[i]+=v; m2[i]+=(v*v);} } END { printf("%i %i %g %g %g %i", m1[0], ell, hopping, epsilon, sigma, dim);  for (i=1; i<=8; i++) { m1[i]/=m1[0]; m2[i]/=m1[0]; printf(" %g %g", m1[i], sqrt((m2[i]-m1[i]*m1[i])/(m1[0]-1))); } printf("\n");} '
done


