WDIR=`pwd`

for i in 108 144 192 256 320 400 500
do

cd $WDIR/$i/nve/

sbatch run.slurm
# cp $WDIR/108/nve/nve.lammps $WDIR/$i/nve/

done