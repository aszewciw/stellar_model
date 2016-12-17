#!/usr/bin/bash

rm ../data/part_a.dat
time ./bin/part_a

python plot_results.py