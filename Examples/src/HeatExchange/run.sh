#!/bin/bash

make
printf "Executing the program with L2 norm\n\n"
./main -p parameters1.pot
printf "\n\n\n\n\n"
printf "Executing the program with H1 norm\n\n"
./main -p parameters2.pot
printf "\n\n\n\n\n"
# gnuplot myplot.gnu
