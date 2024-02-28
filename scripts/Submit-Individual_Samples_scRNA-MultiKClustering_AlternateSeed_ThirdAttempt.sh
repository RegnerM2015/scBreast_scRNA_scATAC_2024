#!/usr/bin/env bash
#SBATCH --job-name Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 4g
#SBATCH --partition allnodes
#SBATCH --output scripts/Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt.out 

# Store sample/container IDs
declare -a StringArray=("4B146L")

# Iterate the string array using for loop, and submit R script for each sample/container id
for i in ${StringArray[@]}
do
    jobfile=scripts/Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt-${i}.sh
    cat <<EOF > $jobfile
#!/usr/bin/env bash
#SBATCH --job-name Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt-${i}
#SBATCH -c 32
#SBATCH --mem 96g
#SBATCH --partition allnodes
#SBATCH --output scripts/Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt-${i}.out 

# Activate conda env that has singularity installed
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate singularity

# Run R script in docker container 
singularity exec --bind /datastore \\
                 --home $PWD \\
                 docker://regnerm/scbreast_2023:1.0.5 \\
                 R CMD BATCH '--args ${i}' scripts/Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt.R scripts/Individual_Samples_scRNA-MultiKClustering_AlternateSeed_ThirdAttempt-${i}.Rout
EOF
    chmod +x ${jobfile}
    sbatch --export=ALL ${jobfile}
done
