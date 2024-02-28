#!/usr/bin/env bash
#SBATCH --job-name Submit-Individual_CellLine_Samples_scRNA-SCSubtype_Classification
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 4g
#SBATCH --partition allnodes
#SBATCH --output scripts/Submit-Individual_CellLine_Samples_scRNA-SCSubtype_Classification.out 

# Store sample/container IDs
declare -a StringArray=("MCF7 T47D HCC1143 SUM149PT")

# Iterate the string array using for loop, and submit R script for each sample/container id
for i in ${StringArray[@]}
do
    jobfile=scripts/Individual_CellLine_Samples_scRNA-SCSubtype_Classification-${i}.sh
    cat <<EOF > $jobfile
#!/usr/bin/env bash
#SBATCH --job-name Individual_CellLine_Samples_scRNA-SCSubtype_Classification-${i}
#SBATCH -c 32
#SBATCH --mem 112g
#SBATCH --partition allnodes
#SBATCH --output scripts/Individual_CellLine_Samples_scRNA-SCSubtype_Classification-${i}.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore \\
                 --home $PWD \\
                 docker://regnerm/scbreast_2023:1.0.5 \\
                 R CMD BATCH '--args ${i}' scripts/Individual_CellLine_Samples_scRNA-SCSubtype_Classification.R scripts/Individual_CellLine_Samples_scRNA-SCSubtype_Classification-${i}.Rout
EOF
    chmod +x ${jobfile}
    sbatch --export=ALL ${jobfile}
done
