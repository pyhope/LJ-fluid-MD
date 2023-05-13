#!/bin/bash

if [ -z $1 ]; then
    echo Missing parameters!
    exit 1
fi

mkdir nve
mkdump.py -i nvt.dump -o select.dump -t $1
mv select.dump nve/
cp nve.lammps nve/
cd nve
atomsk select.dump nve.lmp
sed -i '12c 1   39.94800200             # Ar' nve.lmp
cat ../run.slurm | sed -e 's/01:30:00/24:00:00/g' -e 's/nvt.lammps/nve.lammps/g' > ./run.slurm
sbatch run.slurm