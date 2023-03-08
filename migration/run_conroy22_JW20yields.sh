#!/bin/bash
DTD_LIST=(
	"powerlaw" 
    "powerlaw" 
	"plateau" 
	"plateau" 
	"exponential" 
	"exponential" 
	"prompt" 
    "triple" 
)
DTD_PARAMS=(
	"slope=-1.1" 
	"slope=-1.4" 
	"width=0.3_slope=-1.1" 
	"width=1_slope=-1.1" 
	"timescale=1.5" 
	"timescale=3" 
	"peak=0.05_stdev=0.015_timescale=3"
    "" 
)
DTD_NAMES=(
	"powerlaw_slope11" 
	"powerlaw_slope14" 
	"plateau_width300_slope11" 
	"plateau_width1000_slope11" 
	"exponential_timescale15" 
	"exponential_timescale30" 
	"prompt_peak050_stdev015_timescale30" 
    "triple_delay040" 
)
EVOL="conroy22"
MIGR_LIST=("post-process" "diffusion")

for i in ${!DTD_LIST[@]}; do
	for MIGR in ${MIGR_LIST[@]}; do
		DTD=${DTD_LIST[$i]}
		DTD_PARAM=${DTD_PARAMS[$i]}
		NAME=../src/data/migration/$MIGR/$EVOL"_JW20yields"/${DTD_NAMES[$i]}
		echo $NAME
		python simulations.py -f --nstars=8 --migration=$MIGR --evolution=$EVOL --yields="JW20" --RIa=$DTD --RIa-params=$DTD_PARAM --minimum-delay=0.04 --name=$NAME
	done
done
