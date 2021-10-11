#PBS -l walltime=400:00:00
#PBS -l mem=8gb
#PBS -m ae
#PBS -N SNAKEMASTER
#PBS -o logs/workflow.o
#PBS -e logs/workflow.e

cd ${PBS_O_WORKDIR}

module load bbc/snakemake/snakemake-5.28.0

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --use-envmodules -n > logs/runs/workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/runs/workflow_${TIME}.png

snakemake \
--use-envmodules \
--jobs 20 \
--cluster "qsub \
-V \
-l nodes=1:ppn={threads} \
-l mem={resources.mem_gb}gb \
-l walltime={resources.walltime} \
-o logs/runs/ \
-e logs/runs/"


