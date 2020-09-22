DRIVER_LOC=$(cat ../locations/MDI_metadynamics)
LAMMPS_LOC=$(cat ../locations/LAMMPS)

#remove old files
if [ -d work ]; then 
  rm -r work
fi

#create work directory
cp -r data work
cd work


export OMP_NUM_THREADS=1

#launch LAMMPS
${LAMMPS_LOC} -mdi "-role ENGINE -name LAMMPS -method TCP -port 6022 -hostname localhost" -in in.spcfw_nacl > lammps.out &

#launch driver
${DRIVER_LOC} -mdi "-role DRIVER -name driver -method TCP -port 6022" &

wait
