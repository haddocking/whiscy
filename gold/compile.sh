#!/bin/bash
CFLAGS="DEBUG=0"

g++ -D ${CFLAGS} -o whiscy whiscy.cpp pamcalc.cpp -O3
