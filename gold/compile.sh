#!/bin/bash
CFLAGS="DEBUG=0"

g++ -D ${CFLAGS} -o whiscy whiscy.cpp pamcalc.cpp -O3
g++ -D ${CFLAGS} -o parasmooth parasmooth.cpp -O3
g++ -D ${CFLAGS} -o consadjust consadjust.cpp -O3
