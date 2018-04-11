#!/bin/bash

ulimit -s  102400

#gcc -c ./inp.c -I ./
#ar crv inp.a inp.o

g++ -c ./VIPLimiter.cpp -g
g++ -c ./SphericalmovingGRP.cpp -g
g++ -o SphericalmovingGRP.out ./SphericalmovingGRP.o ./VIPLimiter.o -lm

#RUN

./SphericalmovingGRP.out

exit 0
