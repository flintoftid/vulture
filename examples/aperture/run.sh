#!/bin/sh

runjob()
{

  base="$1"
  xmin="$2"
  emax="$3"

  # Create picture of mesh.
  gvulture -p "${base}.mesh"
  gnuplot mesh.gnp
  mv mesh.eps "${base}-mesh.eps"

  # Run solver and save ASCII observer data.
  vulture "${base}.mesh"
  mv eh_inc1_td.asc "${base}-eh_inc1_td.asc"
  mv eh_inc2_fd.asc "${base}-eh_inc2_fd.asc"
  mv eh_trans1_td.asc "${base}-eh_trans1_td.asc"
  mv eh_trans2_fd.asc "${base}-eh_trans2_fd.asc"

  # Create data files for animation sequence.
  cat <<EOF > process.dat
CE  Vulture example: Aperture in infinite plane
1
 0 60 1 0 40 1 20 20 1
1 200
0 5.99586e+11 2.99793e+09
0.000999998
 0 60 0 40 20 20 3
EOF

  mxplane

  # Render frames to bitmaps.
  for frame in `seq 1 200`
  do

    num=`printf "%06d" ${frame}`

    echo "set ticslevel 0" > frame.gnp
    echo "set style data lines" >> frame.gnp  
    echo "set param" >> frame.gnp
    echo "set zrange [-${emax}:${emax}]" >> frame.gnp
    echo "set xlabel 'x (cells)'" >> frame.gnp
    echo "set ylabel 'y (cells)'" >> frame.gnp
    echo "set xrange [${xmin}:60]" >> frame.gnp
    echo "set zlabel 'E_z (V/m)' offset 7,5" >> frame.gnp
    echo "set title 'time-step: ${frame}'" >> frame.gnp
    echo "se terminal png font 'Arial,12' enhanced large"  >> frame.gnp
    echo "se output 'frame${num}.png'" >> frame.gnp
    echo "" >> frame.gnp
    echo "splot 'tds${frame}.dat' ti '' ls 1" >> frame.gnp
    echo "" >> frame.gnp

    gnuplot frame.gnp

  done

  # Create pictures of selected frames.
  for frame in 50 100 110 120 130 140 150
  do

    num=`printf "%06d" ${frame}`

    echo "set ticslevel 0" > frame.gnp
    echo "set style data lines" >> frame.gnp  
    echo "set param" >> frame.gnp
    echo "set zrange [-${emax}:${emax}]" >> frame.gnp
    echo "set xlabel 'x (cells)'" >> frame.gnp
    echo "set ylabel 'y (cells)'" >> frame.gnp
    echo "set xrange [${xmin}:60]" >> frame.gnp
    echo "set zlabel 'E_z (V/m)' offset 7,5" >> frame.gnp
    echo "set title 'time-step: ${frame}'" >> frame.gnp
    echo "se terminal post eps enhanced color 'Helvetica' 16"  >> frame.gnp
    echo "se output 'frame${num}.eps'" >> frame.gnp
    echo "" >> frame.gnp
    echo "splot 'tds${frame}.dat' ti '' ls 1" >> frame.gnp
    echo "" >> frame.gnp

    gnuplot frame.gnp

    cp "frame${num}.eps" "${base}-frame${num}.eps"

  done

  # Creatre animation.
  ffmpeg -i 'frame%06d.png' -target pal-dvd "${base}.mpg"

  # Tidy up.
  rm -f frame*.png tds*.dat eh*.asc *.log *.dat frame.gnp plot mesh.gnp

}

# Empty space with plane-wave.
runjob "aperture-empty"         0 1.10

# Solid PEC plane with indicent plane-wave.
runjob "aperture-solid"         0 1.10

# Solid PEC plane with incident and reflected plane-waves.
runjob "aperture-solidreflect"  0 1.10

# PEC plane with aperture.
runjob "aperture"              30 0.01

# Solid PEC plane and equivalanet dipole model of aperture.
runjob "aperture-equivdipole"  30 0.01

# Create comparison plots for tangential magnetic field at aperture.
gnuplot plotfd_Hsc.gnp

# Create comparison plots for transmitted electric field.
gnuplot plottd_Erad.gnp

