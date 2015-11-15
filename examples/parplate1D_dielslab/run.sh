#!/bin/sh

# Plot mesh.
gvulture -p parplate1D_dielslab.mesh
gnuplot mesh.gnp
mv mesh.eps parplate1D_dielslab_mesh.eps

# Run solver.
vulture -v parplate1D_dielslab.mesh

# Extract times response along mesh at all time-steps.
cat <<EOF > process.dat
CE  Vulture example: Parallel plate waveguide with dielectric slab
1
 0 200 1 0 0 1 0 0 1
1 1000
0 5.99585e+10 1.19917e+07
0.01
 0 200 0 0 0 0 2
EOF

mxplane

# Plot spatial response at each time-step to bitmap.
for frame in `seq 1 1000`
do

  num=`printf "%06d" ${frame}`

  echo "se style data lines" > frame.gnp  
  echo "set xrange [0:200]" >> frame.gnp
  echo "set yrange [-1:1]" >> frame.gnp
  echo "set xlabel 'x (cells)'" >> frame.gnp
  echo "set ylabel 'E_y (V/m)'" >> frame.gnp
  echo "se term png font 'Arial,12' enhanced large"  >> frame.gnp
  echo "se ou 'frame${num}.png'" >> frame.gnp
  echo "" >> frame.gnp
  echo "se tit 'time-step: ${frame}'" >> frame.gnp
  echo "se arrow 1 from  30,-1 to  30,1 nohead lt 0" >> frame.gnp 
  echo "se arrow 2 from 170,-1 to 170,1 nohead lt 0" >> frame.gnp 
  echo "se label 1 'SF' at  26,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 2 'TF' at  34,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 3 'TF' at 166,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 4 'SF' at 174,-0.9 left rotate by 90" >> frame.gnp
  echo "se arrow 3 from  50,-1 to  50,1 nohead lt 0" >> frame.gnp 
  echo "se arrow 4 from 150,-1 to 150,1 nohead lt 0" >> frame.gnp 
  echo "se label 5 'Free-space' at  46,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 6 'Dielectric' at  54,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 7 'Dielectric' at 146,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 8 'Free-space' at 154,-0.9 left rotate by 90" >> frame.gnp
  echo "plot 'tds${frame}.dat' us 1:2 ti '' ls 1" >> frame.gnp
  echo "" >> frame.gnp

  gnuplot frame.gnp

done

# Assemble bitmaps into animation.
ffmpeg -i 'frame%06d.png' -target pal-dvd parplate1D_dielslab.mpg

# Plot spatial response every 50 time-steps.
for frame in 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800
do

  num=`printf "%06d" ${frame}`

  echo "se style data lines" > frame.gnp  
  echo "set xrange [0:200]" >> frame.gnp
  echo "set yrange [-1:1]" >> frame.gnp
  echo "set xlabel 'x (cells)'" >> frame.gnp
  echo "set ylabel 'E_y (V/m)'" >> frame.gnp
  echo "se term post eps enhanced color 'Helvetica' 16"  >> frame.gnp
  echo "se ou 'frame${num}.eps'" >> frame.gnp
  echo "" >> frame.gnp
  echo "se tit 'time-step: ${frame}'" >> frame.gnp
  echo "se arrow 1 from  30,-1 to  30,1 nohead lt 0" >> frame.gnp 
  echo "se arrow 2 from 170,-1 to 170,1 nohead lt 0" >> frame.gnp 
  echo "se label 1 'SF' at  26,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 2 'TF' at  34,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 3 'TF' at 166,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 4 'SF' at 174,-0.9 left rotate by 90" >> frame.gnp
  echo "se arrow 3 from  50,-1 to  50,1 nohead lt 0" >> frame.gnp 
  echo "se arrow 4 from 150,-1 to 150,1 nohead lt 0" >> frame.gnp 
  echo "se label 5 'Free-space' at  46,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 6 'Dielectric' at  54,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 7 'Dielectric' at 146,-0.9 left rotate by 90" >> frame.gnp
  echo "se label 8 'Free-space' at 154,-0.9 left rotate by 90" >> frame.gnp
  echo "plot 'tds${frame}.dat' us 1:2 ti '' ls 1" >> frame.gnp
  echo "" >> frame.gnp

  gnuplot frame.gnp

  cp "frame${num}.eps" "parplate1D_dielslab-frame${num}.eps"

done

# Create data file with analytic solution.
octave analytic.m

# Plot scattering parameters.
gnuplot plotfd_Smag.gnp
gnuplot plotfd_Sarg.gnp

# Tidy up.
rm -f frame*.png tds*.dat *.log  frame.gnp plot mesh.gnp *.dat *.asc

exit 0

