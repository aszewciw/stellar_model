#!/usr/bin/bash

rm ../data/output.dat
time ./bin/part_a

python plot_results.py