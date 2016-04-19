#!/bin/bash
for i in {7..9}; do cd-hit -i ../../data/MSU7/all.pep -o cdhit-0.$i.out -c 0.$i -n 5 -d 0 -T 1; done
