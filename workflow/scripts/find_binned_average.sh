#!/usr/bin/env bash
################################################################################
##
## Find average binned methylation for input BED file
##
## Created by:
##   Jacob Morrison
##
## Creation date:
##   Mar 2021
##
## Update notes:
##   Mar 2021 -
##     - Initial creation
##   Nov 2021
##     - replace "filtered.tmp.bed" with hash to prevent problems when running multiple jobs in parallel
##
################################################################################

set -euo pipefail



# Find binned average values
function find_averages {
    tmp=$(echo $RANDOM|md5sum|head -c 20; echo) # get random name for tmp file
    if [[ `file ${infile}` =~ "gzip" ]]; then
        zcat ${infile} |
        awk -v filt=${filter} '{ if ($5 >= filt) { print } }' > $tmp
    else
        awk -v filt=${filter} '{ if ($5 >= filt) { print } }' \
        ${infile} > $tmp
    fi

    # Find average beta value across bins
    bedtools makewindows -w ${window} \
        -g ${reffai} | \
    sort -k1,1 -k2,2n | \
    bedtools map \
        -prec 2 \
        -a - \
        -b $tmp \
        -c 4 -o mean |
    gzip > ${otfile}.gz

    rm -f $tmp
}

# Print helpful usage information
function usage {
    >&2 echo -e "\nUsage: find_binned_average.sh [-h,--help] reference window filter infile outfile"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\treference    : FAI-index for reference to create windows from"
    >&2 echo -e "\twindow       : Window size for binned average finding"
    >&2 echo -e "\tfilter       : Coverage filter to apply before finding averages"
    >&2 echo -e "\tinfile       : BISCUIT CG BED file to calculate methylation average from"
    >&2 echo -e "\toutfile      : Name of BED file to write output to"
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "\t-h,--help    : Print help message and exit\n"
}

# Process command line arguments
OPTS=$(getopt \
    --options h \
    --long help \
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
if [[ $# -ne 5 ]]; then
    >&2 echo "$0: Missing inputs"
    usage
    exit 1
fi

# Fill required positional arguments
reffai="${1}"
window="${2}"
filter="${3}"
infile="${4}"
otfile="${5}"

# Do some checks on the given inputs
if [[ ! -f ${reffai} ]]; then
    >&2 echo "Doesn't appear like you provided an existing FAI file: ${reffai}"
    exit 1
fi

if [[ ! "${reffai}" =~ .fai$ ]]; then
    >&2 echo "Doesn't appear like you provided a viable FAI file: ${reffai}"
    exit 1
fi

if ! [[ "${window}" =~ ^[0-9]+$ ]]; then
    >&2 echo "window needs to be an integer: ${window}"
    exit 1
fi

if ! [[ ${filter} =~ ^[0-9]+$ ]]; then
    >&2 echo "filter needs to be an integer: ${filter}"
    exit 1
fi

if [[ ! -f ${infile} ]]; then
    >&2 echo "Doesn't appear that infile exists: ${infile}"
    exit 1
fi

find_averages
