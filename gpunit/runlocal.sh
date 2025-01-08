#!/bin/bash



docker run -v $PWD:$PWD  genepattern/audit:0.4 /build/audit/run_audit.sh $PWD/skyline_export/skyline-export-data.xls  skyline-export-data-results 0.00001  0.2 FALSE TRUE FALSE


