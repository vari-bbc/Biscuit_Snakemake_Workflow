#PBS -l walltime=400:00:00
#PBS -l mem=8gb
#PBS -m ae
#PBS -N SNAKEMASTER
#PBS -o logs/workflow.o
#PBS -e logs/workflow.e

# Function to pull down config settings from config/config.yaml
# Pulled from StackOverflow:
#     https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script
function parse_yaml {
    local prefix=$2
    local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
    awk -F$fs '{
        indent = length($1)/2;
        vname[indent] = $2;
        for (i in vname) {if (i > indent) {delete vname[i]}}
        if (length($3) > 0) {
            vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
            printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
        }
    }'
}

cd ${PBS_O_WORKDIR}

# Extract output directory from config.yaml, set to "." for default
#     behavior (i.e., writing to same directory as Snakefile)
eval $(parse_yaml config/config.yaml "CONF_")
if [[ ${CONF_output_directory} == "" ]]; then
    : ${CONF_output_directory:="."}
fi
mkdir -p ${CONF_output_directory}/logs/runs

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --use-conda -n > ${CONF_output_directory}/logs/runs/workflow_${TIME}.txt
snakemake --dag | dot -Tpng > ${CONF_output_directory}/logs/runs/workflow_${TIME}.png

# Default to using conda, if using environment modules, then replace --use-conda with --use-envmodules
# Note, this requires downloading mamba (conda install -n base -c conda-forge mamba)
snakemake \
    --use-conda \
    --jobs 20 \
    --cluster "qsub \
    -V \
    -l nodes=1:ppn={threads} \
    -l mem={resources.mem_gb}gb \
    -l walltime={resources.walltime} \
    -o ${CONF_output_directory}/logs/runs/ \
    -e ${CONF_output_directory}/logs/runs/"


