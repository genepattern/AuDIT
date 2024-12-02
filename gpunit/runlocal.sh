#!/bin/bash



docker run -v $PWD:$PWD genepattern/audit:0.1 java -DR_HOME=/usr/local/lib64/R -Dr_flags="--no-save --quiet --slave --no-restore" -cp /build RunR /build/audit/AuDIT.R parse.cmdline /build/audit/ $PWD/skyline_export/skyline-export-data.xls  skyline-export-data-results 0.00001  0.2 FALSE TRUE FALSE


