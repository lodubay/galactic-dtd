#!/bin/bash
for DTD in "powerlaw" "powerlaw_steep" "powerlaw_broken" "exponential" "exponential_long" "bimodal"
do
	for MIGR in "diffusion" "post-process"
	do
		python simulations.py -f --nstars=8 --migration=$MIGR --evolution=insideout --RIa=$DTD --name=outputs/$MIGR/insideout/$DTD
	done
done
for DTD in "powerlaw" "powerlaw_steep" "exponential"
do
	for MIGR in "diffusion" "post-process"
	do
		python simulations.py -f --nstars=8 --migration=$MIGR --evolution=insideout --RIa=$DTD --minimum-delay=0.15 --name=outputs/$MIGR/insideout/$DTD_delayed
	done
done
