#!/bin/bash
CFLAGS=DEBUG

g++ -D ${CFLAGS} -o whiscy whiscy.cpp pamcalc.cpp -O3
