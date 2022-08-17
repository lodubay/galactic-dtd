#!/bin/bash
DTD_LIST=(
	"powerlaw" 
	"powerlaw" 
	"plateau" 
	"plateau" 
	"exponential" 
	"exponential" 
	"prompt"
)
DTD_PARAMS=(
	"slope=-1.1" 
	"slope=-1.4" 
	"width=0.3_slope=-1.1" 
	"width=1_slope=-1.1" 
	"timescale=1.5" 
	"timescale=3" 
	"peak=0.05_stdev=0.015_timescale=3"
)
DTD_NAMES=(
	"powerlaw/slope11" 
	"powerlaw/slope14" 
	"plateau/width300_slope11" 
	"plateau/width1000_slope11" 
	"exponential/timescale15" 
	"exponential/timescale30" 
	"prompt/peak050_timescale30"
)
EVOL_LIST=("insideout" "lateburst" "twoinfall" "conroy22")
MIGR_LIST=("diffusion" "post-process")

for EVOL in ${EVOL_LIST[@]}; do
	for i in ${!DTD_LIST[@]}; do
		for MIGR in ${MIGR_LIST[@]}; do
			DTD=${DTD_LIST[$i]}
			DTD_PARAM=${DTD_PARAMS[$i]}
			NAME=../src/data/migration/$MIGR/$EVOL/${DTD_NAMES[$i]}
			echo $NAME
			python simulations.py -f --nstars=8 --migration=$MIGR --evolution=$EVOL --RIa=$DTD --RIa-params=$DTD_PARAM --minimum-delay=0.04 --name=$NAME
		done
	done
done
