#!/bin/bash

for i in {1..19}
do
    grep chr$i $1 > $2chr"$i"_indels.bed
done

grep chrX $1 > $2chrX_indels.bed
grep chrX $1 > $2chrY_indels.bed
