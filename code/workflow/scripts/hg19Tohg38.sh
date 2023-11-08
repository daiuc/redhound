#!/bin/bash


echo './liftover.sh my_hg19.bed my_hg38.bed my.unmapped'
liftOver $1  ~/software/liftover/hg19ToHg38.over.chain $2 $3
