#!/env/bin/bash

# SBATCH --job-name=hwas
# SBATCH --partition=$partion
# SBATCH --account=$account
# SBATCH --qos=$qos
# SBATCH --nodes=$nodes
# SBATCH --time=$time
# SBATCH --time=$logdir/
#


module purge
module load hwas

