#!/usr/bin/env bash
################################################################################
##
## Create SCGS, SCGW, and WCGW context BED files
##
## Notes:
##     1.) bedtools, awk, and parallel must be in PATH for script to work
##
## Created by:
##     Jacob Morrison (based on code from Wanding Zhou)
##
## Creation date:
##     Feb 2022
##
## Update notes:
##     Feb 2022 -
##         Initial creation
################################################################################

set -euo pipefail

# Check for bedtools, awk, and parallel in PATH
function check_path {
    if [[ `which bedtools 2>&1 > /dev/null` ]]; then
        >&2 echo "bedtools does not exist in PATH"
        exit 1
    else
        >&2 echo "Using bedtools found at: `which bedtools`"
    fi
    if [[ `which awk 2>&1 > /dev/null` ]]; then
        >&2 echo "awk does not exist in PATH"
        exit 1
    else
        if awk --version | grep -q GNU; then
            >&2 echo "Using GNU awk found at: `which awk`"
        else
            >&2 echo "It doesn't appear you are using GNU awk"
            >&2 echo "Try adding GNU awk at the front of PATH"
            exit 1
        fi
    fi
    if [[ `which parallel 2>&1 > /dev/null` ]]; then
        >&2 echo "parallel does not exist in PATH"
        >&2 echo "Make sure to add GNU parallel to PATH"
        exit 1
    else
        if parallel --version | grep -q GNU; then
            >&2 echo "Using GNU parallel found at: `which parallel`"
        else
            >&2 echo "It doesn't appear you are using GNU parallel."
            >&2 echo "Try adding GNU parallel at the front of PATH"
            exit 1
        fi
    fi
}

expand_windows() {
    # Load chromosome sizes to help with keeping extracted ranges within bounds of the chromosomes
    chr_sizes=""
    while read line; do
        chr="$(echo ${line} | awk '{ print $1 }')"
        siz="$(echo ${line} | awk '{ print $2 }')"

        if [ -z "$chr_sizes" ]; then
            chr_sizes="${chr},${siz}"
        else
            chr_sizes="${chr_sizes};${chr},${siz}"
        fi
    done < $1

    awk -v sizes="${chr_sizes}" \
    'BEGIN { 
        OFS="\t"

        # Create a key-value dictionary of chromosome sizes
        split(sizes, tmp1, ";")
        for (i in tmp1) {
            split(tmp1[i], tmp2, ",")
            chr_lengths[tmp2[1]] = tmp2[2]
        }
    } {
        chr = $1
        if ($2 - 35 < 0) { beg = 0 }
        else { beg = $2 - 35 }
        if ($3 + 35 > chr_lengths[chr]) { end = chr_lengths[chr] }
        else { end = $3 + 35 }

        print chr, beg, end
    }' $2
}
export -f expand_windows

find_context() {
    awk \
    'BEGIN { OFS = "\t" } {
        # Extract position
        match($1, /([^:]*):([0-9]*)-([0-9]*)/, position)
        chr = position[1]
        beg = position[2] + 35
        end = position[3] - 35

        # Reference sequence
        a1=toupper($2)
        a2=a1
        a3=a1

        # 4-base sequence and context
        rseq_4 = substr(a1,35,4)
        if (rseq_4~/[CG]CG[CG]/) { context = "SCGS" }
        else if (rseq_4~/[AT]CG[AT]/) { context = "WCGW" }
        else { context = "SCGW" }

        n_cpgs = gsub(/CG/,"",a1)                # number of CpGs in string
        gc_per = gsub(/[CG]/,"",a2) / length(a3) # GC-content of string
        rseq_6 = substr(a3,34,6)                 # 6-base sequence
        rseq_2 = substr(a3,36,2)                 # 2-base sequence

        # Need to have an easy way to extract CpGs into specific context and number of neighboring CpGs files
        if (n_cpgs == 1) { tag = "0" }
        else if (n_cpgs == 2) { tag = "1" }
        else if (n_cpgs == 3) { tag = "2" }
        else { tag = "3" }

        if (!(rseq_4~/N/)) print chr, beg, end, n_cpgs-1, gc_per, rseq_4, rseq_6, context, rseq_2, context"_"tag
    }' $1;
}
export -f find_context

create_files() {
    check_path

    if [ ! -d ${outdir} ]; then
        mkdir -p ${outdir}
    fi

    cd ${outdir}

    # Split up CpGs for quicker processing
    zcat ${cpgbed} | split -d -l 1000000 - cpgsplit_

    # Retrieve sequences around CpGs from genome FASTA file
    export genome
    parallel -j ${thread} 'expand_windows ${genome}.fai {} | bedtools getfasta -bed - -fi ${genome} -tab -fo cpgcontext_{}.tab' ::: cpgsplit_*

    # Annotate CpGs
    parallel -j ${thread} 'find_context {} > {.}.context' ::: cpgcontext_*.tab

    # Extract CpGs into specific context and number of neighboring CpG files
    grep --no-filename WCGW_0 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > wcgw_0_neighbors.bed.gz
    grep --no-filename WCGW_1 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > wcgw_1_neighbors.bed.gz
    grep --no-filename WCGW_2 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > wcgw_2_neighbors.bed.gz
    grep --no-filename WCGW_3 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > wcgw_3p_neighbors.bed.gz
    grep --no-filename SCGW_0 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgw_0_neighbors.bed.gz
    grep --no-filename SCGW_1 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgw_1_neighbors.bed.gz
    grep --no-filename SCGW_2 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgw_2_neighbors.bed.gz
    grep --no-filename SCGW_3 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgw_3p_neighbors.bed.gz
    grep --no-filename SCGS_0 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgs_0_neighbors.bed.gz
    grep --no-filename SCGS_1 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgs_1_neighbors.bed.gz
    grep --no-filename SCGS_2 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgs_2_neighbors.bed.gz
    grep --no-filename SCGS_3 *.context | cut -f1-3,5-7 | sort -k1,1 -k2,2n -k3,3n | gzip > scgs_3p_neighbors.bed.gz

    rm cpgsplit_* cpgcontext_*
}

usage() {
    >&2 echo -e "\nUsage: create_context_beds.sh [-h,--help] [-o,--outdir] genome cpgbed\n"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\tgenome       : Path to reference FASTA file used in creating CpG BED"
    >&2 echo -e "\tcpgbed       : CpG BED file containing CpG locations in genome\n"
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "\t-h,--help    : Print help message and exit"
    >&2 echo -e "\t-o,--outdir  : Output directory [DEFAULT: assets]"
    >&2 echo -e "\t-t,--threads : Number of threads to use [DEFAULT: 1]"
}

# Initialize default variable values
outdir="assets"
thread=1

# Process command line arguments
OPTS=$(getopt \
    --options ho:t: \
    --long help,outdir:,threads: \
    --name "$(basename "$0")" \
    -- "$@"
)
eval set -- ${OPTS}

while true; do
    case "$1" in
        -h|--help )
            usage
            exit 0
            ;;
        -o|--outdir )
            outdir="$2"
            shift 2
            ;;
        -t|--threads )
            thread="$2"
            shift 2
            ;;
        -- )
            shift
            break
            ;;
        * )
            >&2 echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Make sure there are the correct number of inputs
if [[ $# -ne 2 ]]; then
    >&2 echo "$0: Missing inputs"
    usage
    exit 1
fi

# Fill required positional arguments
genome=$1
cpgbed=$2

# Input checks
if [[ ! -f "${genome}.fai" ]]; then
    >&2 echo "Cannot locate fai-indexed reference: ${genome}.fai"
    >&2 echo "Please provide a viable path to the reference genome FASTA file."
    exit 1
fi

if [[ ! -f "${cpgbed}" ]]; then
    >&2 echo "Cannot locate CpG BED file: ${cpgbed}"
    >&2 echo "Please provide an existing CpG BED file"
    exit 1
fi

>&2 echo -e "Inputs:"
>&2 echo -e "\tOutput directory : ${outdir}"
>&2 echo -e "\tNumber of threads: ${thread}"
>&2 echo -e "\tGenome Reference : ${genome}"
>&2 echo -e "\tCpG Location BED : ${cpgbed}"

create_files
