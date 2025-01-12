#!/bin/bash

#load required modules
#module load jre/8.121
#module load R/4.3.1-foss-2022b

project='MetaPept'
queue='samuels'
threads=20
mem_per_thread=1000

bsub -q ${queue} -J ${project} -eo logs/${project}.error -oo logs/${project}.out -n ${threads} -R "rusage[mem=${mem_per_thread}] span[ptile=${threads}]" "bash MetaPept.sh"
