#!/bin/bash

WDIR=`pwd`
for i in 144 192 256 320 400 500
do
cd $WDIR/$i
cp $WDIR/108/in.lammps .
done

sed -i 's/0 3 0 3 0 3/0 4 0 3 0 3/' $WDIR/144/in.lammps
sed -i 's/0 3 0 3 0 3/0 4 0 4 0 3/' $WDIR/192/in.lammps
sed -i 's/0 3 0 3 0 3/0 4 0 4 0 4/' $WDIR/256/in.lammps
sed -i 's/0 3 0 3 0 3/0 5 0 4 0 4/' $WDIR/320/in.lammps
sed -i 's/0 3 0 3 0 3/0 5 0 5 0 4/' $WDIR/400/in.lammps
sed -i 's/0 3 0 3 0 3/0 5 0 5 0 5/' $WDIR/500/in.lammps
