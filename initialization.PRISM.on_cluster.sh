module load jre/8.121
bsub -J MetaPept -q samuels -eo logs/MetaPept.ini.error -oo logs/MetaPept.ini.out -R "rusage[mem=20000]" "bash initialization.PRISM.sh"
