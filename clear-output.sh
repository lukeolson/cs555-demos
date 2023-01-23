#!/bin/bash

for f in *.ipynb
do
jupyter nbconvert --clear-output --inplace $f
done
