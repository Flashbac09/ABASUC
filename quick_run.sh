#!/bin/bash

g++ -O3 source/*.cpp -o abasuc -L/home/jzw2/mylib/OpenBLAS-0.3.23 -lopenblas;
./abasuc;
