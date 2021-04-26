#!/usr/bin/env bash
#Jordan Hughey
#Liu Group

#this script runs the HiC snp significant contact pipeline

#this is a script that is ran over PBS and runs multiple scripts

#First blocks to get opts from command line
ARGUMENT_LIST=(
    "resolution"
    "resolution_short"
    "tissue"
    "window"
    "cutoff"
)

# read arguments
opts=$(getopt \
    --longoptions "$(printf "%s:," "${ARGUMENT_LIST[@]}")" \
    --name "$(basename "$0")" \
    --options "" \
    -- "$@"
)

eval set --$opts

while [[ $# -gt 0 ]]; do
    case "$1" in
        --resolution)
            resolution=$2
            shift 2
            ;;

        --resolution_short)
            resolution_short=$2
            shift 2
            ;;

        --tissue)
            tissue=$2
            shift 2
            ;;

        --window)
            window=$2
            shift 2
            ;;

        --cutoff)
            cutoff=$2
            shift 2
            ;;

        *)
            break
            ;;
    esac
done



#input for get opts
res=${resolution}
res_short=${resolution_short}
tiss=${tissue}
window=${window}
cutoff=${cutoff}

wd=`pwd`
 
qsub -v res=${res},res_short=${res_short},tiss=${tiss},window=${window},cutoff=${cutoff} -d ${wd} -N Sig_act_${tiss} ./scripts/activity_snp_pipeline.pbs
#qsub -v res=${r},tissue=${t},window=${w},cutoff=${c} -d ${wd} -N Sig_act_${t} activity_snp_pipeline_v2.pbs

