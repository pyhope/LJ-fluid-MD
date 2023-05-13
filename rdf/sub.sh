WDIR=`pwd`

for i in 108 144 192 256 320 400 500
do

cd $WDIR/$i/
sbatch run.slurm

done