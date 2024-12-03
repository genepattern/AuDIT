#!/bin/bash

export R_HOME=/usr/local/lib64/R
#  docker run --rm -v $PWD:$PWD -w $PWD genepattern/audit:0.2 java -DR_HOME=<R_HOME> \"-Dr_flags=--no-save --quiet --slave --no-restore\" -cp /build RunR 

docker run --rm -v $PWD:$PWD -w $PWD genepattern/audit:0.2 java -Dr_flags="--no-save --quiet --slave --no-restore"  -DR_HOME=/usr/local/lib64/R -cp /build RunR /build/audit/AuDIT.R parse.cmdline /build/audit/ $PWD/gpunit/skyline_export/skyline-export-data.xls results .00001 0.2 FALSE TRUE FALSE


