#!/bin/bash
echo 'run_genome_analysis.py' $1
tmux new -d -s "run_genome_analysis_$1" ./runRGA.sh &> "run_genome_analysis_$1.log"
