#!/bin/bash
DTD_LIST = ("powerlaw" "powerlaw" "plateau" "plateau" "exponential" "exponential" "prompt")
DTD_PARAMS = ("{slope=-1.1}" "{slope=-1.4}" "{width=0.3, slope=-1.1}" "{width=1, slope=-1.1}" "{timescale=1.5}" "{timescale=3}" "{peak=0.05, stdev=0.015, timescale=3}")
EVOL_LIST = ("insideout" "lateburst")
MIGR_LIST = ("diffusion" "post-process")
TAU_STAR = "conroy22"

for EVOL in ${EVOL_LIST[@]}; do
	for i in ${!DTD_LIST[@]}; do
		for MIGR in ${MIGR_LIST[@]}; do
			DTD = ${DTD_LIST[$i]}
			DTD_PARAM = ${DTD_PARAMS[$i]}
			NAME = ../src/data/migration/$MIGR/$EVOL/$TAU_STAR/$DTD/$DTD_PARAM
			echo $NAME
			EVOL_FULL = $EVOL"_ifrmode"
			python simulations.py -f --nstars=8 --migration=$MIGR --evolution=$EVOL_FULL --RIa=$DTD --RIa-kwargs=$DTD_PARAM --minimum-delay=0.04 --name=$NAME