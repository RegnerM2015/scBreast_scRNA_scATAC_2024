#!/usr/bin/env bash
#SBATCH --job-name Submit-Individual_Samples_scRNA-inferCNV_CancerCellDetection
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 4g
#SBATCH --partition allnodes
#SBATCH --output scripts/Submit-Individual_Samples_scRNA-inferCNV_CancerCellDetection.out 

# Store sample/container IDs
declare -a StringArray=("35A4AL 35EE8L 3821AL 3B3E9L 3C7D1L 3D388L 3FCDEL 43E7BL 43E7CL 44F0AL 45CB0L 4C2E5L 4D0D2L")

# Iterate the string array using for loop, and submit R script for each sample/container id
for i in ${StringArray[@]}
do
    jobfile=scripts/Individual_Samples_scRNA-inferCNV_CancerCellDetection-${i}.sh
    cat <<EOF > $jobfile
#!/usr/bin/env bash
#SBATCH --job-name Individual_Samples_scRNA-inferCNV_CancerCellDetection-${i}
#SBATCH -c 32
#SBATCH --mem 64g
#SBATCH --partition allnodes
#SBATCH --output scripts/Individual_Samples_scRNA-inferCNV_CancerCellDetection-${i}.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore \\
                 --home $PWD \\
                 docker://regnerm/scbreast_2023:1.0.5 \\
                 R CMD BATCH '--args ${i}' scripts/Individual_Samples_scRNA-inferCNV_CancerCellDetection.R scripts/Individual_Samples_scRNA-inferCNV_CancerCellDetection-${i}.Rout
EOF
    chmod +x ${jobfile}
    sbatch --export=ALL ${jobfile}
done
