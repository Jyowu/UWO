#!/bin/bash
#SBATCH --account=def-cdegroot
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=1024M
#SBATCH --time=4-10:00


module load openfoam/12

decomposePar

sed -i "s/maxCo\s*[0-9.]*/maxCo           2/" system/controlDict
sed -i "s/writeInterval\s*[0-9.]*/writeInterval   5/" system/controlDict
sed -i "s/endTime\s\+[0-9.]*/endTime         50/" system/controlDict
srun foamRun -parallel

# Segment 3: Co=1, t=8 -> 10
sed -i "s/maxCo\s*[0-9.]*/maxCo           1/" system/controlDict
sed -i "s/writeInterval\s*[0-9.]*/writeInterval   5/" system/controlDict
sed -i "s/endTime\s\+[0-9.]*/endTime         70/" system/controlDict
srun foamRun -parallel

