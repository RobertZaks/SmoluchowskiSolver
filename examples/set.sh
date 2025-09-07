#!/bin/sh
for alpha in 0.001 0.005 0.01 0.03 0.02 0.05
do
    ./smolsolver "$alpha" > /dev/null
    mv functional.csv "res/functional${alpha}.csv" 
done
