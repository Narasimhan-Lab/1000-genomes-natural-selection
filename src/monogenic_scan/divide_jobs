#!/bin/bash

total=1150477
i=0
chunk=$((total/$1))

while(($i <= ($total - $chunk)))
do
	echo $i
	sbatch submit_job $i $chunk
	((i=i+chunk))
done
covered=$((chunk*$1))
remaining=$((total- covered))
sbatch single_job $covered $remaining
