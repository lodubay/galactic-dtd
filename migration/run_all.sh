#!/bin/bash
for DTD in "powerlaw" "powerlaw_steep" "powerlaw_broken" "exponential" "exponential_long" "bimodal"
do
	for EVOL in "insideout" "lateburst" "insideout_conroy22" "lateburst_conroy22"
	do
		for MIGR in "diffusion" "post-process"
		do
			NAME=../src/data/migration/$MIGR/$EVOL/$DTD
			echo $NAME
			python simulations.py -f --nstars=8 --migration=$MIGR --evolution=$EVOL --RIa=$DTD --name=$NAME
		done
	done
done
for DTD in "powerlaw" "powerlaw_steep" "exponential"
do
	for EVOL in "insideout" "lateburst" "insideout_conroy22" "lateburst_conroy22"
	do
		for MIGR in "diffusion" "post-process"
		do
			NAME=../src/data/migration/$MIGR/$EVOL/$DTD\_delayed
			echo $NAME
			python simulations.py -f --nstars=8 --migration=$MIGR --evolution=$EVOL --RIa=$DTD --minimum-delay=0.15 --name=$NAME
		done
	done
done
