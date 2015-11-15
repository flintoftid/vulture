#!/bin/sh

for i in aperture enclosure_se free_space parallel_plate parplate1D_dielslab planewave
do

  cd "${i}"
  ./run.sh
  cd ..

done

