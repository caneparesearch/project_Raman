#!/bin/bash

echo -e "\n\n" >> temporary 
sed -n '/FINAL OPTIMIZED GEOMETRY - DIMENSIONALITY OF THE SYSTEM      3/,/T = ATOM BELONGING TO THE ASYMMETRIC UNIT/p' *.out >> temporary
#cat temporary
sed '1, 2d; $d' temporary
rm temporary
