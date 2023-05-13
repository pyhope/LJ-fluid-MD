#!/bin/bash

sed -i 's/0 3 0 3 0 3/0 4 0 3 0 3/' ./144/nvt.lammps
sed -i 's/0 3 0 3 0 3/0 4 0 4 0 3/' ./192/nvt.lammps
sed -i 's/0 3 0 3 0 3/0 4 0 4 0 4/' ./256/nvt.lammps
sed -i 's/0 3 0 3 0 3/0 5 0 4 0 4/' ./320/nvt.lammps
sed -i 's/0 3 0 3 0 3/0 5 0 5 0 4/' ./400/nvt.lammps
sed -i 's/0 3 0 3 0 3/0 5 0 5 0 5/' ./500/nvt.lammps
