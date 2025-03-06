#!/bin/bash

#SBATCH --job-name=test
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH -c 4
#SBATCH --mem-per-cpu=2G
#SBATCH --gres=disk:1000
#SBATCH --time=24:00:00

export MYEXECUTABLE=$SAGA_HOME/build/getKinBinnedAsym
export OUTDIR=$SAGA_HOME/tutorials/results_1D
export YAML=$OUTDIR/args.yaml

echo $MYEXECUTABLE
echo $OUTDIR
echo $YAML

cd $OUTDIR
ls -lrth
pwd
$MYEXECUTABLE $YAML
echo DONE
