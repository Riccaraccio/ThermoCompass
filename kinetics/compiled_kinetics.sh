#!/bin/bash

echo "Compiling All the Kinetics Modules ..."

for dir in */;
do
  echo "Compiling $dir"
  cd $dir
  OpenSMOKEpp_CHEMKIN_PreProcessor.sh --input ../input.dic > /dev/null
  cd ..
done

echo "All Kinetics Modules Compiled Successfully!"
