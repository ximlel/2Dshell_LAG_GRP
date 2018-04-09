#!/bin/bash

ulimit -s  102400

#gcc -c ./inp.c -I ./
#ar crv inp.a inp.o

g++ -c ./SphericalmovingGRP.cpp -g
g++ -o SphericalmovingGRP.out ./SphericalmovingGRP.o -lm

#RUN

./SphericalmovingGRP.out

exit 0
