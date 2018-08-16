#!/bin/bash

rm -rf protdist *.o
gcc -o protdist protdist.c -O3 -lm
